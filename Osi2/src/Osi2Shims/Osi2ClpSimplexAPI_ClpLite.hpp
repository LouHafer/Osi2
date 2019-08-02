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

#include "Osi2ParamBEAPI.hpp"

namespace Osi2 {

class DynamicLibrary ;

/*! \brief Proof of concept API.

  Capable of loading and solving a problem. Period.
*/
class ClpSimplexAPI_ClpLite : public ClpSimplexAPI {

public:
  /// Constructor
  ClpSimplexAPI_ClpLite(DynamicLibrary *libClp) ;

  /// Destructor
  ~ClpSimplexAPI_ClpLite() ;

  /// Read an mps file from the given filename
  int readMps(const char *filename, bool keepNames = false,
	      bool ignoreErrors = false) ;

  /// General solve algorithm
  int initialSolve() ;

private:
  /*! \name Dynamic object management information */
  //@{
    /// Dynamic library handle
    DynamicLibrary *libClp_ ;
    /// Clp object
    Clp_Simplex *clpC_ ;
  //@}

} ;

}  // end namespace Osi2

#endif
