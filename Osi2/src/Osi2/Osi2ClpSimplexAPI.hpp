#ifndef Osi2ClpSimplexAPI_HPP
# define Osi2ClpSimplexAPI_HPP

#include <vector>

#include "Coin_C_defines.h"

namespace Osi2 {

class RunParamsAPI ;

/*! \brief Provide the ClpSimplex `C' interface

  This %API exposes the portion of the ClpSimplex class provided by the Clp
  `C' interface.
*/
class ClpSimplexAPI : public API {

public:

  /// API ID
  inline static const char *getAPIIDString() { return ("ClpSimplex") ; }

  /*! \brief variable status

    First 4 match CoinWarmStartBasis,
    isFixed means fixed at lower bound and out of basis
  */
  enum Status {
    isFree = 0x00,
    basic = 0x01,
    atUpperBound = 0x02,
    atLowerBound = 0x03,
    superBasic = 0x04,
    isFixed = 0x05
  } ;

  /// fake bounds for dual
  enum FakeBound {
    noFake = 0x00,
    lowerFake = 0x01,
    upperFake = 0x02,
    bothFake = 0x03
  } ;

  /*! \name Constructors, destructor, and copy

    Because construction is handled through ControlAPI, we don't need a
    constructor. Arguably, we should not need a destructor, either.

    \todo See if it's possible to get rid of the destructor; at least, make it
    protected.
  */
  //@{

  /// Destructor
  virtual ~ClpSimplexAPI() { }
  //@}

  /*! \name Methods to load, modify, and retrieve a problem */
  //@{
  /*! \brief Load a problem

    The matrix is given in a standard column major ordered format (without
    gaps).
    The constraints on the rows are given by lower and upper bounds. If a
    pointer is 0 then the following values are the default:
    <ul>
      <li> <code>colub</code>: all columns have upper bound infinity
      <li> <code>collb</code>: all columns have lower bound 0
      <li> <code>rowub</code>: all rows have upper bound infinity
      <li> <code>rowlb</code>: all rows have lower bound -infinity
      <li> <code>obj</code>: all variables have 0 objective coefficient
    </ul>
  */
  void loadProblem(const int numcols, const int numrows,
    const CoinBigIndex *start, const int *index,
    const double *value,
    const double *collb, const double *colub,
    const double *obj,
    const double *rowlb, const double *rowub) ;

  /// read quadratic part of the objective (the matrix part)
  void loadQuadraticObjective(const int numberColumns,
      const CoinBigIndex *start, const int *column, const double *element) ;

  /// Get integer information
  char *integerInformation() const ;
  /// Add integer information
  void copyInIntegerInformation(const char *information) ;
  /// Drop integer information
  void deleteIntegerInformation() ;

  /// Resize rim part of model
  void resize(int newNumberRows, int newNumberColumns) ;

  /// Get number of rows
  int getNumRows() ;
  /// Add rows
  void addRows(int number, const double *rowLower, const double *rowUpper,
    const CoinBigIndex *rowStarts, const int *columns,
    const double *elements) ;
  /// Delete rows
  void deleteRows(int number, const int *which) ;

  /// Get number of columns
  int getNumCols() ;
  /// Add columns
  void addColumns(int number,
    const double *columnLower, const double *columnUpper,
    const double *objective,
    const CoinBigIndex *columnStarts, const int *rows,
    const double *elements) ;
  /// Delete columns
  void deleteColumns(int number, const int *which) ;

  /// Get row lower bounds
  double *rowLower() const ;
  /** Change row lower bounds */
  void chgRowLower(const double *rowLower);

  /// Get row upper bounds
  double *rowUpper() const ;
  /** Change row upper bounds */
  void chgRowUpper(const double *rowUpper);

  /// Get column lower bounds
  double *columnLower() const ;
  /** Change column lower bounds */
  void chgColumnLower(const double *columnLower);

  /// Get column upper bounds
  double *columnUpper() const ;
  /** Change column upper bounds */
  void chgColumnUpper(const double *columnUpper);
  
  /// Get objective coefficients
  double *objective() const ;
  /** Change objective coefficients */
  void chgObjCoefficients(const double *objIn);

  /// Number of elements in matrix
  CoinBigIndex getNumElements() const ;
  /// Column starts in matrix
  const CoinBigIndex *getVectorStarts() const ;
  /// Column lengths in matrix
  const int *getVectorLengths() const ;
  /// Row indices in matrix
  const int *getIndices() const ;
  /// element values in matrix
  const double *getElements() const ;
  /// Modify one element of a matrix
  void modifyCoefficient(int row, int column, double newElement,
    bool keepZero = false) ;

  /// Check if status array exists
  bool statusExists() const ;
  /// Get address of status array (char[numberRows+numberColumns])
  unsigned char *statusArray() const ;
  /// Copy in a status vector
  void copyinStatus(const unsigned char *statusArray) ;
  /// Get status for a variable
  ClpSimplexAPI::Status getColumnStatus(int sequence) ;
  /*! \brief Set status for a variable 

    Also force the value if the status is at bound.
  */
  void setColumnStatus(int sequence, ClpSimplexAPI::Status newStatus) ;
  /// Get status for a row
  int getRowStatus(int sequence);
  /*! \brief Set status for a row 

    Also force the value if the status is at bound.
  */
  void setRowStatus(int sequence, ClpSimplexAPI::Status newStatus) ;

  /// Copy in row & column names
  void copyNames(const std::vector<std::string> &rowNames,
    const std::vector<std::string> &columnNames) ;
  /// Drop names
  void dropNames() ;
  /// length of names (0 means no names available)
  int lengthNames() ;
  /// Get row name (size of name must be at least #lengthNames()+1)
  void rowName(int iRow,char *name) ;
  /// Set row name (name should be null-terminated)
  void setRowName(int iRow,char *name) ;
  /// Get column name (size of name must be at least #lengthNames()+1)
  void columnName(int iCol,char *name) ;
  /// Set column name (name should be null-terminated)
  void setColumnName(int iCol,char *name) ;

  //@}

  /*! \name Methods to read or write a problem file */
  //@{
  /// Read an MPS file from the given filename
  virtual int readMps(const char *filename,
    bool keepNames = false,
    bool ignoreErrors = false) = 0 ; 

  /*! \brief Write an MPS file to the given filename

    Format type is 0 = normal, 1 = extra or 2 = hex.
    Number across is 1 or 2.
    Use objSense = -1D to flip the objective function around.
  */
  int writeMps(const char *filename, int formatType, int numberAcross,
    double objSense) ;

  /// Save the model to a file
  int saveModel(const char *fileName) ;
  /// Restore model from a file
  int restoreModel(const char *fileName) ;
  //@}

  /*! \name Parameter Gets & Sets */
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
  virtual int numberIterations() const = 0 ;
  virtual void setNumberIterations(int value) = 0 ;
  /// Maximum number of iterations
  int maximumIterations() ;
  void setMaximumIterations(int value) ;
  /// Maximum time in seconds (timed from call to setMaximumSeconds)
  int maximumSeconds() ;
  void setMaximumSeconds(double value) ;

  /*! \brief Direction of optimisation

    -1.0  maximise
     0.0  ignore
     1.0  minimise
  */
  virtual double objSense() const = 0 ;
  virtual void setObjSense(double value) = 0 ;

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

  /// Get solve algorithm
  int algorithm() ;
  /// Set solve algorithm
  void setAlgorithm(int value) ;

  /*! \brief Get small element value
    Elements less that this are set to zero. Default 1.0e-20.
  */
  double getSmallElementValue() ;
  /// Set small element value
  void setSmallElementValue(double value) ;
  //@}

  /*! \name Bulk parameter manipulation

    These methods work with a RunParamsAPI object to maintain a complete set
    of parameters that can be applied in one operation.
  */
  //@{
  /// Load a runParams object with parameters.
  virtual void exposeParams(RunParamsAPI &runParams) const = 0 ;

  /// Set ClpSimplex parameters from a RunParamsAPI object.
  virtual void loadParams(RunParamsAPI &runParams) = 0 ;
  //@}

  /*! \name Information about the solution */
  //@{
  /*! \brief Status of problem:

       0 - optimal
       1 - primal infeasible
       2 - dual infeasible
       3 - stopped on iterations etc
       4 - stopped due to errors
  */
  int status() ;
  /** Set problem status */
  void _setProblemStatus(int problemStatus) ;
  /*! \brief Secondary status of problem - may get extended

    0 - none
    1 - primal infeasible because dual limit reached
    2 - scaled problem optimal - unscaled has primal infeasibilities
    3 - scaled problem optimal - unscaled has dual infeasibilities
    4 - scaled problem optimal - unscaled has dual and primal infeasibilities
  */
  int secondaryStatus() ;
  void setSecondaryStatus(int status) ;
  /// Problem is primal feasible
  bool primalFeasible() const ;
  /// Problem is dual feasible
  bool dualFeasible() const ;
  /// Are there numerical difficulties?
  bool isAbandoned() ;
  /// Return true if hit maximum iterations or time
  bool hitMaximumIterations() ;
  /// Iteration limit reached?
  bool isIterationLimitReached() ;
  /// Is optimality proven?
  bool isProvenOptimal() ;
  /// Is primal infeasibility proven?
  bool isProvenPrimalInfeasible() ;
  /// Is dual infeasibility proven?
  bool isProvenDualInfeasible() ;
  /// Is the given primal objective limit reached

  /// Objective function value
  double objectiveValue() ;
  /// Primal row solution
  double primalRowSolution() ;
  /// Primal column solution
  double primalColumnSolution() ;
  /// Set the primal column solution
  void setPrimalColumnSolution(double *input) ;

  /*! \brief Unbounded ray (null if bounded)

    Use #freeRay to free the returned ray.
  */
  double *unboundedRay() const ;
  /// Dual row solution
  double dualRowSolution() ;
  /// Dual column solution
  double dualColumnSolution() ;
  /*! Infeasible ray (null if feasible)

    Use #freeRay to free the returned ray.
  */
  double *infeasibilityRay(bool fullRay = false) const ;
  /*! \brief Free a ray

    Use this method to free rays returned by #unboundedRay and
    #infeasibilityRay.
  */
  void freeRay(double *ray) ;

  /// Number of primal infeasibilities
  int numberPrimalInfeasibilities() const ;
  /// Sum of primal infeasibilities
  double sumPrimalInfeasibilities() const ;

  /// Number of dual infeasibilities
  int numberDualInfeasibilities() const ;
  /// Sum of dual infeasibilities
  double sumDualInfeasibilities() const ;

  /*! Check the solution

    Sets sum of infeasibilities, etc.
  */
  void checkSolution() ;
  //@}
  
  /*! \name User pointer
  
    Reserved for use by the application, for whatever purpose.
  */
  //@{
  /// Set the user pointer
  void setUserPointer(void *pointer) ;
  /// Get the user pointer
  void *getUserPointer() ;
  //@}

  /*! \name Message handling

    Callbacks are handled by ONE function.
  */
  //@{
  /// Pass in callback function
  // void registerCallBack(ClpSimplexAPI::Callback userCallBack) ;
  /// Unset callback function
  void clearCallBack() ;
  /// Set the log level
  void setLogLevel(int value) ;
  /// Get the log level
  int logLevel() ;
  //@}


  /*! \name Functions most useful to user */
  //@{
  /// Default initial solve
  virtual int initialSolve() = 0 ;
  /*! \brief General solve algorithm which can do presolve.

     See  ClpSolve.hpp for options
  */
  // int initialSolve(ClpSolve &options) ;
  /// Dual initial solve
  int initialDualSolve() ;
  /// Primal initial solve
  int initialPrimalSolve() ;
  /// Barrier initial solve
  int initialBarrierSolve() ;
  /// Barrier initial solve, not to be followed by crossover
  int initialBarrierNoCrossSolve() ;
  /*! \brief Dual algorithm

    See ClpSimplexDual.hpp for method.
    ifValuesPass==2 just does values pass and then stops.
  */
  int dual(int ifValuesPass = 0) ;
  /*! \brief Primal algorithm

   See ClpSimplexPrimal.hpp for method.
   ifValuesPass==2 just does values pass and then stops.
  */
  int primal(int ifValuesPass = 0) ;

  /*! \brief Crash a basis

    At present, just aimed at dual simplex. Return codes are
    -2 - dual preferred, crash basis created
    -1 - dual preferred, slack basis preferred
     0 - basis going in was not all slack
     1 - primal preferred, slack basis preferred
     2 - primal preferred, crash basis created

   If the gap between bounds <= \p gap, variables can be flipped.

   The values for pivot are
     0 - no pivoting allowed, choice of algorithm only
     1 - simple pivoting
     2 - mini iterations
  */
  int crash(double gap,int pivot) ;
  //@}

} ;

}    // end namespace Osi2

#endif

