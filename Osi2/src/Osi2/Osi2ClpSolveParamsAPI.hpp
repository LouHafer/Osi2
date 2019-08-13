#ifndef Osi2ClpSolveParamsAPI_HPP
# define Osi2ClpSolveParamsAPI_HPP

/*! \file Osi2ClpSolveParamsAPI.hpp

  Provides an %API corresponding to the ClpSolve class.
*/

#include <vector>

#include "Coin_C_defines.h"

namespace Osi2 {

/*! \brief Provide the ClpSolveParams interface

  This %API exposes an object that can hold solve time parameters used to
  control the clp solver. It combines the functionality of the ClpSolve class
  with that of the individual parameter set/get methods from ClpSimplex.
*/
class ClpSolveParamsAPI : public API {

public:

  /// API ID
  inline static const char *getAPIIDString() { return ("ClpSolveParams") ; }

  /*! \name Constructors, destructor, and copy */
  //@{
  /// Destructor
  virtual ~ClpSolveParamsAPI() { }
  //@}

  /*! \name Parameter Gets & Sets

    Methods formerly provided by the ClpSimplex class.
  */
  //@{
  /// Primal zero tolerance
  virtual double primalTolerance() const = 0 ;
  virtual void setPrimalTolerance(double value) = 0 ;

  /// Dual zero tolerance
  virtual double dualTolerance() const = 0 ;
  virtual void setDualTolerance(double value) = 0 ;

  /// Dual objective limit
  double dualObjectiveLimit() ;
  void setDualObjectiveLimit(double value) ;

  /*! \brief Select scaling mode

    0 - off
    1 - equilibrium
    2 - geometric
    3 - auto
    4 - dynamic (implemented?)
  */
  void scaling(int mode) ;
  /// Get scaling mode
  int scalingFlag() ;

  /// Constant offset for objective
  double objectiveOffset() ;
  void setObjectiveOffset(double value) ;

  /// Get problem name (copied to array)
  int problemName(int maxNumberCharacters, char *array) ;
  /// Set problem name to array (null-terminated string)
  int setProblemName(int maxNumberCharacters, char *array) ;

  /// Number of iterations
  int numberIterations() ;
  void setNumberIterations(int value) ;
  /// Maximum number of iterations
  int maximumIterations() ;
  void setMaximumIterations(int value) ;
  /// Maximum time in seconds (timed from call to setMaximumSeconds)
  int maximumSeconds() ;
  void setMaximumSeconds(double value) ;

  /*! \brief Direction of optimisation

    -1  maximise
     0  ignore
     1  minimise
  */
  double getObjSense() ;
  void setObjSense(double value) ;

  /// Dual bound
  double dualBound() ;
  void setDualBound(double value) ;
  
  /// Infeasibility cost
  double infeasibilityCost() ;
  void setInfeasibilityCost(double value) ;

  /*! \brief Problem perturbation

     50 - switch on perturbation
    100 - auto perturb if takes too long (1.0e-6 largest nonzero)
    101 - problem is perturbed
    102 - do not try perturbing again

    Default is 100.
  */
  int perturbation() ;
  void setPerturbation(int value) ;

  /// Get algorithm
  int algorithm() ;
  /// Set algorithm
  void setAlgorithm(int value) ;

  /*! \brief Get small element value
    Elements less that this are set to zero. Default 1.0e-20.
  */
  double getSmallElementValue() ;
  /// Set small element value
  void setSmallElementValue(double value) ;
  //@}

} ;

}    // end namespace Osi2

#endif

