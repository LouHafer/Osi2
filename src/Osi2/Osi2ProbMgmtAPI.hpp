/*
  Copyright 2011 Lou Hafer, Matt Saltzman
  This code is licensed under the terms of the Eclipse Public License (EPL)

  $Id$
*/

#ifndef Osi2ProbMgmtAPI_HPP
#define Osi2ProbMgmtAPI_HPP

#include "Osi2API.hpp"

/*! \brief Proof of concept API.

*/

namespace Osi2 {

class ProbMgmtAPI : public API {

public:

  /// Return the "ident" string for the ProbMgmt %API.
  inline static const char *getAPIIDString () { return ("ProbMgmtAPI") ; }

  /// Virtual destructor
  virtual ~ProbMgmtAPI() {}

  /// Read an mps file from the given filename
  virtual
  int readMps(const char *filename, bool keepNames = false,
              bool ignoreErrors = false) = 0 ;

  /// Solve an lp
  virtual int initialSolve() = 0 ;

  /// API capability
  inline int getAPIs(const char **&idents) {
    apiName_ = getAPIIDString() ;
    idents = &apiName_ ;
    return (1) ;
  }

private:

  /// Intermediate storage for returning the API name
  const char *apiName_ ;

} ;

}  // end namespace Osi2

#endif

