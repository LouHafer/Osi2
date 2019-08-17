/*
  Copyright 2011 Lou Hafer, Matt Saltzman
  This code is licensed under the terms of the Eclipse Public License (EPL)

  $Id$
*/
/*! \file Osi2ClpSimplexAPI_ClpLite.hpp
    \brief The clp `lite' implementation of Osi2::ClpSimplexAPI
*/
#ifndef Osi2ClpSimplexAPI_ClpLite_HPP
#define Osi2ClpSimplexAPI_ClpLite_HPP

#include "Osi2API.hpp"
#define COIN_EXTERN_C
#include "Osi2ClpSimplexAPI.hpp"

#include "Osi2ParamBEAPI_Imp.hpp"

namespace Osi2 {

class DynamicLibrary ;

/*! \brief Proof of concept API.

  Capable of loading and solving a problem. Period.
*/
class ClpSimplexAPI_ClpLite : public ClpSimplexAPI {

public:
/*! \name Constructors, destructor, and copy */
//@{
  /// Constructor
  ClpSimplexAPI_ClpLite(DynamicLibrary *libClp) ;

  /// Destructor
  ~ClpSimplexAPI_ClpLite() ;
//@}

/*! \name Methods to read or write a problem file */
//@{
  /// Read an mps file from the given filename
  int readMps(const char *filename, bool keepNames = false,
	      bool ignoreErrors = false) ;
//@}

/*! \name Parameter gets & sets */
//@{
  /// Primal zero tolerance
  double primalTolerance() const ;
  void setPrimalTolerance(double val) ;

  /// Dual zero tolerance
  double dualTolerance() const ;
  void setDualTolerance(double val) ;

  /// Number of iterations
  int numberIterations() const ;
  void setNumberIterations(int val) ;

  /*! \brief objective sense

    -1 to maximise, 1 to minimise
  */
  double objSense() const ;
  void setObjSense(double val) ;
//@}

/*! \name Bulk parameter manipulation

  Methods for getting and setting parameters in a way that matches OSI2
  standard parameter management. Also methods to work with a RunParamsAPI
  object to get and set many parameters in a single operation.
*/
//@{
  /// Get / set an integer parameter
  int getIntParam(std::string name) const ;
  void setIntParam(std::string name, int val) ;
  /// Get / set a double parameter
  double getDblParam(std::string name) const ;
  void setDblParam(std::string name, double val) ;
  /// Get / set a string parameter
  std::string getStrParam(std::string name) const ;
  void setStrParam(std::string name, std::string val) ;

  /// Fill a RunParamsAPI object with parameters and current values
  void exposeParams(RunParamsAPI &runParams) const ;
  /// Load parameters from a RunParamsAPI object
  void loadParams(RunParamsAPI &runParams) ;
//@}

/*! \name Functions must useful to user */
//@{
  /// General solve algorithm
  int initialSolve() ;
//@}

/*! \name API management & enquiry */
//@{
  inline int getAPIs(const char **&idents)
  { return (paramMgr_.getAPIs(idents)) ; }

  inline void *getAPIPtr (const char *ident)
  { return (paramMgr_.getAPIPtr(ident)) ; }

private:
  /*! \name Dynamic object management information */
  //@{
    /// Parameter management object
    ParamBEAPI_Imp<ClpSimplexAPI_ClpLite> paramMgr_ ;
    /// Dynamic library handle
    DynamicLibrary *libClp_ ;
    /// Clp object
    Clp_Simplex *clpC_ ;
  //@}

} ;

}  // end namespace Osi2

#endif

