/*
  Copyright 2011 Lou Hafer, Matt Saltzman
  This code is licensed under the terms of the Eclipse Public License (EPL)

  $Id$
*/
/*! \file Osi2ClpSimplexAPI_ClpLite.cpp
    \brief Method definitions for Osi2ClpSimplexAPI_ClpLite

  Method definitions for ClpSimplexAPI_ClpLite, an implementation of the
  ClpSimplex API using the `lite' clp shim that dynamically loads methods
  on demand.
*/

#include <iostream>

#include "Osi2Config.h"
#include "Osi2nullptr.hpp"

#include "Osi2DynamicLibrary.hpp"
#include "Osi2RunParamsAPI.hpp"

#include "Osi2ClpSimplexAPI_ClpLite.hpp"

namespace {

// Templates for common function signatures

#include "Osi2CFuncTemplates.cpp"

/*
  Static maps for know ClpSimplex integer, double, and string parameters. The
  first string is the name exposed through Osi2 standard parameter interfaces,
  the second the name of the C interface get method, the third the name of
  the C interface set method.
*/

KnownMemberMap integerSimplexParams = {
  mapEntry3("iter count","Clp_numberIterations","Clp_setNumberIterations"),
  mapEntry3("iter max","Clp_maximumIterations","Clp_setMaximumIterations"),
  mapEntry3("scaling","Clp_scalingFlag","Clp_scaling"),
  mapEntry3("perturbation","Clp_perturbation","Clp_setPerturbation"),
  mapEntry3("algorithm","Clp_algorithm","Clp_setAlgorithm")
} ;
  

KnownMemberMap doubleSimplexParams = {
  mapEntry3("primal zero tol","Clp_primalTolerance","Clp_setPrimalTolerance"),
  mapEntry3("dual zero tol","Clp_dualTolerance","Clp_setDualTolerance"),
  mapEntry3("dual obj lim",
  	    "Clp_dualObjectiveLimit","Clp_setDualObjectiveLimit"),
  mapEntry3("obj offset","Clp_objectiveOffset","Clp_setObjectiveOffset"),
  mapEntry3("obj sense",
  	    "Clp_optimizationDirection","Clp_setOptimizationDirection"),
  mapEntry3("infeas cost","Clp_infeasibilityCost","Clp_setInfeasibilityCost"),
  mapEntry3("small elem",
  	    "Clp_getSmallElementValue","Clp_setSmallElementValue"),
  mapEntry3("dual bound","Clp_dualBound","Clp_setDualBound")
} ;


KnownMemberMap stringSimplexParams = {
  mapEntry3("problem name","Clp_problemName","Clp_setProblemName")
} ;


/*
  Static map for known presolve transform names. The first string is the name
  exposed through the parameter interface, the second the name of the
  corresponding Clp_C_Interface ClpSolve method.
*/

typedef std::map<std::string,std::string>KnownXformMap ;
KnownXformMap knownXforms = {
  {std::string("return if infeasible"),std::string("infeasibleReturn")},
  {std::string("dual"),std::string("doDual")},
  {std::string("singleton"),std::string("doSingleton")},
  {std::string("doubleton"),std::string("doDoubleton")},
  {std::string("tripleton"),std::string("doTripleton")},
  {std::string("tighten"),std::string("doTighten")},
  {std::string("forcing"),std::string("doForcing")},
  {std::string("implied free"),std::string("doImpliedFree")},
  {std::string("dupcol"),std::string("doDupcol")},
  {std::string("duprow"),std::string("doDuprow")},
  {std::string("singleton column"),std::string("doSingletonColumn")},
  {std::string("kill small"),std::string("doKillSmall")}
} ;

}    // end anonymous namespace


namespace Osi2 {

typedef ClpSimplexAPI_ClpLite CSA_CL ;

/*
  Constructor
*/
CSA_CL::ClpSimplexAPI_ClpLite (DynamicLibrary *libClp)
    : paramMgr_(ParamBEAPI_Imp<CSA_CL>(this)),
      libClp_(libClp),
      clpC_(nullptr)
{
  std::string errStr ;
/*
  Look up Clp_newModel, coerce it into a suitable type, and invoke it to get a
  new Clp_Simplex object.
*/
  typedef Clp_Simplex *(*ClpFactory)() ;
  ClpFactory factory =
      libClp_->getFunc<ClpFactory>(std::string("Clp_newModel"),errStr) ;
  if (factory != nullptr) {
    clpC_ = factory() ;
  }
/*
  Register that we are a ClpSimplex API and that we support parameter
  management.
*/
  paramMgr_.addAPIID(ClpSimplexAPI::getAPIIDString(),this) ;
  paramMgr_.addAPIID(ParamBEAPI::getAPIIDString(),&paramMgr_) ;
/*
  Register the parameters exposed through parameter management
*/
  typedef ParamBEAPI_Imp<CSA_CL>::ParamEntry_Imp<double> DPE ;
  typedef ParamBEAPI_Imp<CSA_CL>::ParamEntry_Imp<int> IPE ;
  ParamBEAPI_Imp<CSA_CL>::ParamEntry *paramEntry ;
  paramEntry = new DPE("primal tolerance",
  		       &CSA_CL::primalTolerance,
  		       &CSA_CL::setPrimalTolerance) ;
  paramMgr_.addParam("primal tolerance",paramEntry) ;
  paramEntry = new DPE("dual tolerance",
  		       &CSA_CL::dualTolerance,
  		       &CSA_CL::setDualTolerance) ;
  paramMgr_.addParam("dual tolerance",paramEntry) ;
  paramEntry = new DPE("obj sense",
  		       &CSA_CL::objSense,
  		       &CSA_CL::setObjSense) ;
  paramMgr_.addParam("dual tolerance",paramEntry) ;
  paramEntry = new IPE("iter count",
  		       &CSA_CL::numberIterations,
  		       &CSA_CL::setNumberIterations) ;
  paramMgr_.addParam("iter count",paramEntry) ;
  return ;
}

CSA_CL::~ClpSimplexAPI_ClpLite ()
{
  std::string errStr ;
/*
  Look up Clp_deleteModel, coerce it into a suitable type, and invoke it to
  destroy the Clp_Simplex object.
*/
  typedef void (*ClpDeleteModelFunc)(Clp_Simplex *clp) ;
  ClpDeleteModelFunc deleteModel =
     libClp_->getFunc<ClpDeleteModelFunc>("Clp_deleteModel",errStr) ;
  if (deleteModel != nullptr) {
    deleteModel(clpC_) ;
    clpC_ = nullptr ;
  }
}

/*
  Read a problem file in mps format.
*/
int CSA_CL::readMps (const char *filename, bool keepNames,
                              bool ignoreErrors)
{
  std::string errStr ;

  typedef int (*ClpReadMpsFunc)(Clp_Simplex *,const char *,int,int) ;
  static ClpReadMpsFunc readMps = nullptr ;

  if (readMps == nullptr) {
    readMps = libClp_->getFunc<ClpReadMpsFunc>("Clp_readMps",errStr) ;
  }
  int retval = -1 ;
  if (readMps != nullptr) {
    retval = readMps(clpC_,filename,keepNames,ignoreErrors) ;
    if (retval) {
	std::cout
	    << "Failure to read " << filename << ", error " << retval
	    << "." << std::endl ;
    } else {
	std::cout
	    << "Read " << filename << " without error." << std::endl ;
    }
  }
  return (retval) ;
}


/*
  Solve a problem
*/
int CSA_CL::initialSolve ()
{
  std::string errStr ;

  typedef int (*ClpInitSolFunc)(Clp_Simplex *) ;
  static ClpInitSolFunc initialSolve = nullptr ;

  if (initialSolve == nullptr) {
    initialSolve = libClp_->getFunc<ClpInitSolFunc>("Clp_initialSolve",errStr) ;
  }
  int retval = -1 ;
  if (initialSolve != nullptr) {
      retval = initialSolve(clpC_) ;
      if (retval < 0) {
	  std::cout
	      << "Solve failed; error " << retval << "." << std::endl ;
      } else {
	  std::cout
	      << "Solved; return status " << retval << "." << std::endl ;
      }
  }
  return (retval) ;
}


/* Individual parameter set / get methods */

double CSA_CL::primalTolerance () const
{
  return (simpleGetter<Clp_Simplex,double>
	      (libClp_,clpC_,"Clp_primalTolerance")) ;
}
void CSA_CL::setPrimalTolerance (double val)
{
  simpleSetter<Clp_Simplex,double>
      (libClp_,clpC_,"Clp_setPrimalTolerance",val) ;
}

double CSA_CL::dualTolerance () const
{
  return (simpleGetter<Clp_Simplex,double>
	      (libClp_,clpC_,"Clp_dualTolerance")) ;
}
void CSA_CL::setDualTolerance (double val)
{
  simpleSetter<Clp_Simplex,double>
      (libClp_,clpC_,"Clp_setDualTolerance",val) ;
}

double CSA_CL::objSense () const
{
  return (simpleGetter<Clp_Simplex,double>
	      (libClp_,clpC_,"Clp_optimizationDirection")) ;
}
void CSA_CL::setObjSense (double val)
{
  simpleSetter<Clp_Simplex,double>
      (libClp_,clpC_,"Clp_setOptimizationDirection",val) ;
}

int CSA_CL::numberIterations () const
{
  return (simpleGetter<Clp_Simplex,int>
	      (libClp_,clpC_,"Clp_numberIterations")) ;
}
void CSA_CL::setNumberIterations (int val)
{
  simpleSetter<Clp_Simplex,int>
      (libClp_,clpC_,"Clp_setNumberIterations",val) ;
}

/*
  Generic parameter set / get methods
*/

int CSA_CL::getIntParam (std::string name) const
{
  return (getMember<Clp_Simplex,int>
	      (integerSimplexParams,name,libClp_,clpC_)) ;
}
double CSA_CL::getDblParam (std::string name) const
{
  return (getMember<Clp_Simplex,double>
	      (doubleSimplexParams,name,libClp_,clpC_)) ;
}
std::string CSA_CL::getStrParam (std::string name) const
{
  return (getMember<Clp_Simplex,std::string>
	      (stringSimplexParams,name,libClp_,clpC_)) ;
}

void CSA_CL::setIntParam (std::string name, int val)
{
  setMember<Clp_Simplex,int>
      (integerSimplexParams,name,libClp_,clpC_,val) ;
}
void CSA_CL::setDblParam (std::string name, double val)
{
  setMember<Clp_Simplex,double>
      (doubleSimplexParams,name,libClp_,clpC_,val) ;
}
void CSA_CL::setStrParam (std::string name, std::string val)
{
  setMember<Clp_Simplex,std::string>
      (stringSimplexParams,name,libClp_,clpC_,val) ;
}

/*
  Methods to work with RunParamsAPI objects.
*/
void CSA_CL::exposeParams (RunParamsAPI &runParams) const
{
  for (KnownMemberMap::const_iterator iter = integerSimplexParams.begin() ;
       iter != integerSimplexParams.end() ;
       iter++) {
    std::string name = iter->first ;
    std::string getter = iter->second.getter_ ;
    int val = simpleGetter<Clp_Simplex,int>(libClp_,clpC_,name) ;
    runParams.addIntParam(name,val) ;
  }
  for (KnownMemberMap::const_iterator iter = doubleSimplexParams.begin() ;
       iter != doubleSimplexParams.end() ;
       iter++) {
    std::string name = iter->first ;
    std::string getter = iter->second.getter_ ;
    double val = simpleGetter<Clp_Simplex,double>(libClp_,clpC_,name) ;
    runParams.addDblParam(name,val) ;
  }
  for (KnownMemberMap::const_iterator iter = stringSimplexParams.begin() ;
       iter != stringSimplexParams.end() ;
       iter++) {
    std::string name = iter->first ;
    std::string getter = iter->second.getter_ ;
    std::string val =
        simpleGetter<Clp_Simplex,std::string>(libClp_,clpC_,name) ;
    runParams.addStrParam(name,val) ;
  }
}

void CSA_CL::loadParams (RunParamsAPI &runParams)
{
  std::vector<std::string> paramNames = runParams.getIntParamIds() ;
  for (std::vector<std::string>::const_iterator iter = paramNames.begin() ;
       iter != paramNames.end() ;
       iter++) {
    std::string name = *iter ;
    int val = runParams.getIntParam(name) ;
    std::string setFunc = integerSimplexParams[name].setter_ ;
    simpleSetter<Clp_Simplex,int>(libClp_,clpC_,setFunc,val) ;
  }
  paramNames = runParams.getDblParamIds() ;
  for (std::vector<std::string>::const_iterator iter = paramNames.begin() ;
       iter != paramNames.end() ;
       iter++) {
    std::string name = *iter ;
    double val = runParams.getDblParam(name) ;
    std::string setFunc = doubleSimplexParams[name].setter_ ;
    simpleSetter<Clp_Simplex,double>(libClp_,clpC_,setFunc,val) ;
  }
  paramNames = runParams.getStrParamIds() ;
  for (std::vector<std::string>::const_iterator iter = paramNames.begin() ;
       iter != paramNames.end() ;
       iter++) {
    std::string name = *iter ;
    std::string val = runParams.getStrParam(name) ;
    std::string setFunc = stringSimplexParams[name].setter_ ;
    simpleSetter<Clp_Simplex,std::string>(libClp_,clpC_,setFunc,val) ;
  }
}

}    // end namespace Osi2

