/*
  Copyright 2011 Lou Hafer, Matt Saltzman
  This code is licensed under the terms of the Eclipse Public License (EPL)

  $Id$
*/

/*! \file Osi2Osi1API.hpp

  This file defines the class Osi2::Osi1API, an API that matches the original
  OsiSolverInterface API.
*/

#ifndef Osi2Osi1API_HPP
#define Osi2Osi1API_HPP

#include <cstdlib>
#include <string>
#include <vector>

#include "CoinMessageHandler.hpp"
#include "CoinPackedVectorBase.hpp"
#include "CoinTypes.h"

#include "OsiCollections.hpp"
#include "OsiSolverParameters.hpp"

#include "Osi2API.hpp"

class CoinPackedMatrix;
class CoinWarmStart;
class CoinSnapshot;
class CoinLpIO;
class CoinMpsIO;

class OsiCuts;
class OsiAuxInfo;
class OsiRowCut;
class OsiRowCutDebugger;
class CoinSet;
class CoinBuild;
class CoinModel;
class OsiSolverBranch;
class OsiSolverResult;
class OsiObject;

namespace Osi2 {


/*! \brief Original OSI interface virtual base class

  This API matches the original OsiSolverInterface API.
*/

class Osi1API : public API  {

public:

  /// Return the "ident" string for the OSI1 %API.
  inline static const char *getAPIIDString () { return ("Osi1API") ; }

  /// API capability
  inline int getAPIs(const char **&idents) {
    apiName_ = getAPIIDString() ;
    idents = &apiName_ ;
    return (1) ;
  }

  /*! \brief Nested class for returning status from the applyCuts method.

    Note that it's impossible for an outside entity to change any of the
    values after construction. The class is const by design.
  */
  class ApplyCutsReturnCode {
  public:
    /// \name Constructors and destructors
    //@{
      /// Default constructor
      ApplyCutsReturnCode():
	 intInconsistent_(0),
	 extInconsistent_(0),
	 infeasible_(0),
	 ineffective_(0),
	 applied_(0) {} 
      /// Constructor with initial values
      ApplyCutsReturnCode(int intIn, int extIn,
      			  int infeas, int ineff, int applied)
       : intInconsistent_(intIn),
	 extInconsistent_(extIn),
	 infeasible_(infeas),
	 ineffective_(ineff),
	 applied_(applied) {} 
      /// Copy constructor
      ApplyCutsReturnCode(const ApplyCutsReturnCode & rhs):
	 intInconsistent_(rhs.intInconsistent_),
	 extInconsistent_(rhs.extInconsistent_),
	 infeasible_(rhs.infeasible_),
	 ineffective_(rhs.ineffective_),
	 applied_(rhs.applied_) {}
      /// Assignment operator
      ApplyCutsReturnCode & operator=(const ApplyCutsReturnCode& rhs)
      { 
	if (this != &rhs) { 
	  intInconsistent_ = rhs.intInconsistent_;
	  extInconsistent_ = rhs.extInconsistent_;
	  infeasible_      = rhs.infeasible_;
	  ineffective_     = rhs.ineffective_;
	  applied_         = rhs.applied_;
	}
	return *this;
      }
      /// Destructor
      virtual ~ApplyCutsReturnCode(){}
    //@}

    /*! \name Accessing return code attributes */
    //@{
      /// Number of logically inconsistent cuts
      inline int getNumInconsistent() const
      {return intInconsistent_;}
      /// Number of cuts inconsistent with the current model
      inline int getNumInconsistentWrtIntegerModel() const
      {return extInconsistent_;}
      /// Number of cuts that cause obvious infeasibility
      inline int getNumInfeasible() const
      {return infeasible_;}
      /// Number of redundant or ineffective cuts
      inline int getNumIneffective() const
      {return ineffective_;}
      /// Number of cuts applied
      inline int getNumApplied() const
      {return applied_;}
    //@}

  private: 
    /*! \name Private methods */
    //@{
      /// Increment logically inconsistent cut counter 
      inline void incrementInternallyInconsistent(){intInconsistent_++;}
      /// Increment model-inconsistent counter
      inline void incrementExternallyInconsistent(){extInconsistent_++;}
      /// Increment infeasible cut counter
      inline void incrementInfeasible(){infeasible_++;}
      /// Increment ineffective cut counter
      inline void incrementIneffective(){ineffective_++;}
      /// Increment applied cut counter
      inline void incrementApplied(){applied_++;}
    //@}

    /// \name Private member data
    //@{
      /// Counter for logically inconsistent cuts
      int intInconsistent_;
      /// Counter for model-inconsistent cuts
      int extInconsistent_;
      /// Counter for infeasible cuts
      int infeasible_;
      /// Counter for ineffective cuts
      int ineffective_;
      /// Counter for applied cuts
      int applied_;
    //@}
  };


  //---------------------------------------------------------------------------

  /// \name Constructors and destructors
  //@{
    
    /** Clone (virtual copy)

      The result of calling clone(false) is defined to be equivalent to
      calling the default constructor Osi1API().
    */
    virtual Osi1API * clone(bool copyData = true) const = 0;
  
    /// Virtual Destructor 
    virtual ~Osi1API () {} ;

    /** Reset the solver interface.

    A call to reset() returns the solver interface to the same state as
    it would have if it had just been constructed by calling the default
    constructor Osi1API() (if there were such a thing).
    */
    virtual void reset() = 0 ;
  //@}

  /// \name Solve methods 
  //@{
    /// Solve initial LP relaxation 
    virtual void initialSolve() = 0; 

    /*! \brief Resolve an LP relaxation after problem modification

      Note the `re-' in `resolve'. initialSolve() should be used to solve the
      problem for the first time.
    */
    virtual void resolve() = 0;
  //@}

  /*! \name Parameter set/get methods

     The set methods return true if the parameter was set to the given value,
     false otherwise. When a set method returns false, the original value (if
     any) should be unchanged.  There can be various reasons for failure: the
     given parameter is not applicable for the solver (e.g., refactorization
     frequency for the volume algorithm), the parameter is not yet
     implemented for the solver or simply the value of the parameter is out
     of the range the solver accepts. If a parameter setting call returns
     false check the details of your solver.

     The get methods return true if the given parameter is applicable for the
     solver and is implemented. In this case the value of the parameter is
     returned in the second argument. Otherwise they return false.

     \note
     There is a default implementation of the set/get
     methods, namely to store/retrieve the given value using an array in the
     base class. A specific solver implementation can use this feature, for
     example, to store parameters that should be used later on. Implementors
     of a solver interface should overload these functions to provide the
     proper interface to and accurately reflect the capabilities of a
     specific solver.

     The format for hints is slightly different in that a boolean specifies
     the sense of the hint and an enum specifies the strength of the hint.
     Hints should be initialised when a solver is instantiated.
     (See OsiSolverParameters.hpp for defined hint parameters and strength.)
     When specifying the sense of the hint, a value of true means to work with
     the hint, false to work against it.  For example,
     <ul>
       <li> \code setHintParam(OsiDoScale,true,OsiHintTry) \endcode
	    is a mild suggestion to the solver to scale the constraint
	    system.
       <li> \code setHintParam(OsiDoScale,false,OsiForceDo) \endcode
	    tells the solver to disable scaling, or throw an exception if
	    it cannot comply.
     </ul>
     As another example, a solver interface could use the value and strength
     of the \c OsiDoReducePrint hint to adjust the amount of information
     printed by the interface and/or solver.  The extent to which a solver
     obeys hints is left to the solver.  The value and strength returned by
     \c getHintParam will match the most recent call to \c setHintParam,
     and will not necessarily reflect the solver's ability to comply with the
     hint.  If the hint strength is \c OsiForceDo, the solver is required to
     throw an exception if it cannot perform the specified action.

     \note
     As with the other set/get methods, there is a default implementation
     which maintains arrays in the base class for hint sense and strength.
     The default implementation does not store the \c otherInformation
     pointer, and always throws an exception for strength \c OsiForceDo.
     Implementors of a solver interface should override these functions to
     provide the proper interface to and accurately reflect the capabilities
     of a specific solver.
  */
  //@{
    //! Set an integer parameter
    virtual bool setIntParam(OsiIntParam key, int value) = 0 ;
    //! Set a double parameter
    virtual bool setDblParam(OsiDblParam key, double value) = 0 ;
    //! Set a string parameter
    virtual bool setStrParam(OsiStrParam key, const std::string & value) = 0 ;
    /*! \brief Set a hint parameter

      The \c otherInformation parameter can be used to pass in an arbitrary
      block of information which is interpreted by the OSI and the underlying
      solver.  Users are cautioned that this hook is solver-specific.

      Implementors:
      The default implementation completely ignores \c otherInformation and
      always throws an exception for OsiForceDo. This is almost certainly not
      the behaviour you want; you really should override this method.
    */
    virtual bool setHintParam(OsiHintParam key, bool yesNo=true,
			      OsiHintStrength strength=OsiHintTry,
			      void *otherInformation = nullptr) = 0 ;
    //! Get an integer parameter
    virtual bool getIntParam(OsiIntParam key, int& value) const = 0 ;
    //! Get a double parameter
    virtual bool getDblParam(OsiDblParam key, double& value) const = 0 ;
    //! Get a string parameter
    virtual bool getStrParam(OsiStrParam key, std::string& value) const = 0 ;
    /*! \brief Get a hint parameter (all information)

	Return all available information for the hint: sense, strength,
	and any extra information associated with the hint.

	Implementors: The default implementation will always set
	\c otherInformation to nullptr. This is almost certainly not the
	behaviour you want; you really should override this method.
    */
    virtual bool getHintParam(OsiHintParam key, bool& yesNo,
			      OsiHintStrength& strength,
			      void *& otherInformation) const = 0 ;
    /*! \brief Get a hint parameter (sense and strength only)

	Return only the sense and strength of the hint.
    */
    virtual bool getHintParam(OsiHintParam key, bool& yesNo,
			      OsiHintStrength& strength) const = 0 ;
    /*! \brief Get a hint parameter (sense only)

	Return only the sense (true/false) of the hint.
    */
    virtual bool getHintParam(OsiHintParam key, bool& yesNo) const = 0 ;
    /*! \brief Copy all parameters in this section from one solver to another

      Note that the current implementation also copies the appData block,
      message handler, and rowCutDebugger. Arguably these should have
      independent copy methods.
    */
    virtual void copyParameters(Osi1API & rhs) = 0 ;

    /** \brief Return the integrality tolerance of the underlying solver.
    
	We should be able to get an integrality tolerance, but
        until that time just use the primal tolerance

	\todo
	This method should be replaced; it's architecturally wrong.  This
	should be an honest dblParam with a keyword. Underlying solvers
	that do not support integer variables should return false for set and
	get on this parameter.  Underlying solvers that support integrality
	should add this to the parameters they support, using whatever
	tolerance is appropriate.  -lh, 091021-
    */
    virtual double getIntegerTolerance() const = 0 ;
  //@}

  /// \name Methods returning info on how the solution process terminated
  //@{
    /// Are there numerical difficulties?
    virtual bool isAbandoned() const = 0 ;
    /// Is optimality proven?
    virtual bool isProvenOptimal() const = 0 ;
    /// Is primal infeasibility proven?
    virtual bool isProvenPrimalInfeasible() const = 0 ;
    /// Is dual infeasibility proven?
    virtual bool isProvenDualInfeasible() const = 0 ;
    /// Is the given primal objective limit reached?
    virtual bool isPrimalObjectiveLimitReached() const = 0 ;
    /// Is the given dual objective limit reached?
    virtual bool isDualObjectiveLimitReached() const = 0 ;
    /// Iteration limit reached?
    virtual bool isIterationLimitReached() const = 0 ;
  //@}

  /*! \name Warm start methods

    Note that the warm start methods return a generic CoinWarmStart object.
    The precise characteristics of this object are solver-dependent. Clients
    who wish to maintain a maximum degree of solver independence should take
    care to avoid unnecessary assumptions about the properties of a warm start
    object.
  */
  //@{
    /*! \brief Get an empty warm start object
      
      This routine returns an empty warm start object. Its purpose is
      to provide a way for a client to acquire a warm start object of the
      appropriate type for the solver, which can then be resized and modified
      as desired.
    */

    virtual CoinWarmStart *getEmptyWarmStart () const = 0 ;

    /** \brief Get warm start information.

      Return warm start information for the current state of the solver
      interface. If there is no valid warm start information, an empty warm
      start object wil be returned.
    */
    virtual CoinWarmStart* getWarmStart() const = 0;
    /** \brief Get warm start information.

      Return warm start information for the current state of the solver
      interface. If there is no valid warm start information, an empty warm
      start object wil be returned.  This does not necessarily create an 
      object - may just point to one.  mustDelete is set true if user
      should delete returned object.
    */
    virtual CoinWarmStart* getPointerToWarmStart(bool & mustDelete) = 0 ;

    /** \brief Set warm start information.
    
      Return true or false depending on whether the warm start information was
      accepted or not.
      By definition, a call to setWarmStart with a null parameter should
      cause the solver interface to refresh its warm start information
      from the underlying solver.
   */
    virtual bool setWarmStart(const CoinWarmStart* warmstart) = 0 ;
  //@}

  /*! \name Hot start methods
  
     Primarily used in strong branching. The user can create a hot start
     object --- a snapshot of the optimization process --- then reoptimize
     over and over again, starting from the same point.

     \note
     <ul>
     <li> Between hot started optimizations only bound changes are allowed.
     <li> The copy constructor and assignment operator should NOT copy any
          hot start information.
     <li> The default implementation simply extracts a warm start object in
          \c markHotStart, resets to the warm start object in
	  \c solveFromHotStart, and deletes the warm start object in
	  \c unmarkHotStart.
	  <em>Actual solver implementations are encouraged to do better.</em>
     </ul>

  */
  //@{
    /// Create a hot start snapshot of the optimization process.
    virtual void markHotStart() = 0 ;
    /// Optimize starting from the hot start snapshot.
    virtual void solveFromHotStart() = 0 ;
    /// Delete the hot start snapshot.
    virtual void unmarkHotStart() = 0 ;
  //@}

  /*! \name Problem query methods

   Querying a problem that has no data associated with it will result in
   zeros for the number of rows and columns, and null pointers from the
   methods that return vectors.

   Const pointers returned from any data-query method are valid as long as
   the data is unchanged and the solver is not called.
  */
  //@{
    /// Get the number of columns
    virtual int getNumCols() const = 0;

    /// Get the number of rows
    virtual int getNumRows() const = 0;

    /// Get the number of nonzero elements
    virtual int getNumElements() const = 0;

    /// Get the number of integer variables
    virtual int getNumIntegers() const = 0 ;

    /// Get a pointer to an array[getNumCols()] of column lower bounds
    virtual const double * getColLower() const = 0;

    /// Get a pointer to an array[getNumCols()] of column upper bounds
    virtual const double * getColUpper() const = 0;

    /*! \brief Get a pointer to an array[getNumRows()] of row constraint senses.

      <ul>
      <li>'L': <= constraint
      <li>'E': =  constraint
      <li>'G': >= constraint
      <li>'R': ranged constraint
      <li>'N': free constraint
      </ul>
    */
    virtual const char * getRowSense() const = 0;

    /*! \brief Get a pointer to an array[getNumRows()] of row right-hand sides

      <ul>
	<li> if getRowSense()[i] == 'L' then
	     getRightHandSide()[i] == getRowUpper()[i]
	<li> if getRowSense()[i] == 'G' then
	     getRightHandSide()[i] == getRowLower()[i]
	<li> if getRowSense()[i] == 'R' then
	     getRightHandSide()[i] == getRowUpper()[i]
	<li> if getRowSense()[i] == 'N' then
	     getRightHandSide()[i] == 0.0
      </ul>
    */
    virtual const double * getRightHandSide() const = 0;

    /*! \brief Get a pointer to an array[getNumRows()] of row ranges.

      <ul>
	  <li> if getRowSense()[i] == 'R' then
		  getRowRange()[i] == getRowUpper()[i] - getRowLower()[i]
	  <li> if getRowSense()[i] != 'R' then
		  getRowRange()[i] is 0.0
	</ul>
    */
    virtual const double * getRowRange() const = 0;

    /// Get a pointer to an array[getNumRows()] of row lower bounds
    virtual const double * getRowLower() const = 0;

    /// Get a pointer to an array[getNumRows()] of row upper bounds
    virtual const double * getRowUpper() const = 0;

    /*! \brief Get a pointer to an array[getNumCols()] of objective
	       function coefficients.
    */
    virtual const double * getObjCoefficients() const = 0;

    /*! \brief Get the objective function sense
    
      -  1 for minimisation (default)
      - -1 for maximisation
    */
    virtual double getObjSense() const = 0;

    /// Return true if the variable is continuous
    virtual bool isContinuous(int colIndex) const = 0;

    /// Return true if the variable is binary
    virtual bool isBinary(int colIndex) const = 0 ;

    /*! \brief Return true if the variable is integer.

      This method returns true if the variable is binary or general integer.
    */
    virtual bool isInteger(int colIndex) const = 0 ;

    /// Return true if the variable is general integer
    virtual bool isIntegerNonBinary(int colIndex) const = 0 ;

    /// Return true if the variable is binary and not fixed
    virtual bool isFreeBinary(int colIndex) const = 0 ; 

    /*! \brief Return an array[getNumCols()] of column types

       - 0 - continuous
       - 1 - binary
       - 2 - general integer
      
      If \p refresh is true, the classification of integer variables as
      binary or general integer will be reevaluated. If the current bounds
      are [0,1], or if the variable is fixed at 0 or 1, it will be classified
      as binary, otherwise it will be classified as general integer.
    */
    virtual const char * getColType(bool refresh=false) const = 0 ;
  
    /// Get a pointer to a row-wise copy of the matrix
    virtual const CoinPackedMatrix * getMatrixByRow() const = 0;

    /// Get a pointer to a column-wise copy of the matrix
    virtual const CoinPackedMatrix * getMatrixByCol() const = 0;

    /*! \brief Get a pointer to a mutable row-wise copy of the matrix.
    
      Returns nullptr if the request is not meaningful (i.e., the OSI will not
      recognise any modifications to the matrix).
    */
    virtual CoinPackedMatrix * getMutableMatrixByRow() const = 0 ;

    /*! \brief Get a pointer to a mutable column-wise copy of the matrix
    
      Returns nullptr if the request is not meaningful (i.e., the OSI will not
      recognise any modifications to the matrix).
    */
    virtual CoinPackedMatrix * getMutableMatrixByCol() const = 0 ;

    /// Get the solver's value for infinity
    virtual double getInfinity() const = 0;
  //@}

  /*! \name Solution query methods */
  //@{
    /// Get a pointer to an array[getNumCols()] of primal variable values
    virtual const double * getColSolution() const = 0;

    /** Get a pointer to an array[getNumCols()] of primal variable values
	guaranteed to be between the column lower and upper bounds.
    */
    virtual const double * getStrictColSolution() = 0 ;

    /// Get pointer to array[getNumRows()] of dual variable values
    virtual const double * getRowPrice() const = 0;

    /// Get a pointer to an array[getNumCols()] of reduced costs
    virtual const double * getReducedCost() const = 0;

    /** Get a pointer to array[getNumRows()] of row activity levels.
    
      The row activity for a row is the left-hand side evaluated at the
      current solution.
    */
    virtual const double * getRowActivity() const = 0;

    /// Get the objective function value.
    virtual double getObjValue() const = 0;

    /** Get the number of iterations it took to solve the problem (whatever
	`iteration' means to the solver).
    */
    virtual int getIterationCount() const = 0;

    /** Get as many dual rays as the solver can provide. In case of proven
	primal infeasibility there should (with high probability) be at least
	one.

	The first getNumRows() ray components will always be associated with
	the row duals (as returned by getRowPrice()). If \c fullRay is true,
	the final getNumCols() entries will correspond to the ray components
	associated with the nonbasic variables. If the full ray is requested
	and the method cannot provide it, it will throw an exception.

	\note
	Implementors of solver interfaces note that the double pointers in
	the vector should point to arrays of length getNumRows() (fullRay =
	false) or (getNumRows()+getNumCols()) (fullRay = true) and they should
	be allocated with new[].

	\note
	Clients of solver interfaces note that it is the client's
	responsibility to free the double pointers in the vector using
	delete[]. Clients are reminded that a problem can be dual and primal
	infeasible.
    */
    virtual std::vector<double*> getDualRays(int maxNumRays,
					     bool fullRay = false) const = 0;

    /** Get as many primal rays as the solver can provide. In case of proven
	dual infeasibility there should (with high probability) be at least
	one.
   
	\note
	Implementors of solver interfaces note that the double pointers in
	the vector should point to arrays of length getNumCols() and they
	should be allocated with new[].

	\note
	Clients of solver interfaces note that it is the client's
	responsibility to free the double pointers in the vector using
	delete[]. Clients are reminded that a problem can be dual and primal
	infeasible.
    */
    virtual std::vector<double*> getPrimalRays(int maxNumRays) const = 0;

    /** Get vector of indices of primal variables which are integer variables 
	but have fractional values in the current solution. */
    virtual OsiVectorInt getFractionalIndices(const double etol=1.e-05) const = 0 ;
  //@}

  /*! \name Methods to modify the objective, bounds, and solution

     For functions which take a set of indices as parameters
     (\c setObjCoeffSet(), \c setColSetBounds(), \c setRowSetBounds(),
     \c setRowSetTypes()), the parameters follow the C++ STL iterator
     convention: \c indexFirst points to the first index in the
     set, and \c indexLast points to a position one past the last index
     in the set.
  
  */
  //@{
    /** Set an objective function coefficient */
    virtual void setObjCoeff( int elementIndex, double elementValue ) = 0;

    /** Set a set of objective function coefficients */
    virtual void setObjCoeffSet(const int* indexFirst,
				const int* indexLast,
				const double* coeffList) = 0 ;

    /** Set the objective coefficients for all columns.

	array [getNumCols()] is an array of values for the objective.
	This defaults to a series of set operations and is here for speed.
    */
    virtual void setObjective(const double * array) = 0 ;

    /** Set the objective function sense.

        Use 1 for minimisation (default), -1 for maximisation.

	\note
	Implementors note that objective function sense is a parameter of
	the OSI, not a property of the problem. Objective sense can be
	set prior to problem load and should not be affected by loading a
	new problem.
    */
    virtual void setObjSense(double s) = 0;
  

    /*! Set a single column lower bound.

	Use -getInfinity() for -infinity.
    */
    virtual void setColLower( int elementIndex, double elementValue ) = 0;
    
    /** Set the lower bounds for all columns.

	\p array[getNumCols()] is an array of values for the lower bounds.
	This defaults to a series of set operations and is here for speed.
    */
    virtual void setColLower(const double * array) = 0 ;

    /** Set a single column upper bound.

	Use getInfinity() for infinity.
    */
    virtual void setColUpper( int elementIndex, double elementValue ) = 0;

    /** Set the upper bounds for all columns.

	array [getNumCols()] is an array of values for the upper bounds.
	This defaults to a series of set operations and is here for speed.
    */
    virtual void setColUpper(const double * array) = 0 ;
    
    
    /** Set a single column lower and upper bound.

	The default implementation just invokes setColLower() and
	setColUpper()
    */
    virtual void setColBounds( int elementIndex,
			       double lower, double upper ) = 0 ;
  
    /** Set the upper and lower bounds of a set of columns.

	The default implementation just invokes setColBounds() over and over
	again.  For each column, boundList must contain both a lower and
	upper bound, in that order.
    */
    virtual void setColSetBounds(const int* indexFirst,
				 const int* indexLast,
				 const double* boundList) = 0 ;

    /** Set a single row lower bound.
	Use -getInfinity() for -infinity. */
    virtual void setRowLower( int elementIndex, double elementValue ) = 0;
    
    /** Set a single row upper bound.
	Use getInfinity() for infinity. */
    virtual void setRowUpper( int elementIndex, double elementValue ) = 0;
  
    /** Set a single row lower and upper bound.
	The default implementation just invokes setRowLower() and
	setRowUpper() */
    virtual void setRowBounds( int elementIndex,
			       double lower, double upper ) = 0 ;

    /** Set the bounds on a set of rows.

	The default implementation just invokes setRowBounds() over and over
	again.  For each row, boundList must contain both a lower and
	upper bound, in that order.
    */
    virtual void setRowSetBounds(const int* indexFirst,
				 const int* indexLast,
				 const double* boundList) = 0 ;
  
  
    /** Set the type of a single row */
    virtual void setRowType(int index, char sense, double rightHandSide,
			    double range) = 0;
  
    /** Set the type of a set of rows.
	The default implementation just invokes setRowType()
	over and over again.
    */
    virtual void setRowSetTypes(const int* indexFirst,
				const int* indexLast,
				const char* senseList,
				const double* rhsList,
				const double* rangeList) = 0 ;

    /** Set the primal solution variable values

	colsol[getNumCols()] is an array of values for the primal variables.
	These values are copied to memory owned by the solver interface
	object or the solver.  They will be returned as the result of
	getColSolution() until changed by another call to setColSolution() or
	by a call to any solver routine.  Whether the solver makes use of the
	solution in any way is solver-dependent.
    */
    virtual void setColSolution(const double *colsol) = 0;

    /** Set dual solution variable values

	rowprice[getNumRows()] is an array of values for the dual variables.
	These values are copied to memory owned by the solver interface
	object or the solver.  They will be returned as the result of
	getRowPrice() until changed by another call to setRowPrice() or by a
	call to any solver routine.  Whether the solver makes use of the
	solution in any way is solver-dependent.
    */
    virtual void setRowPrice(const double * rowprice) = 0;

    /** Fix variables at bound based on reduced cost
    
	For variables currently at bound, fix the variable at bound if the
	reduced cost exceeds the gap. Return the number of variables fixed.

	If justInteger is set to false, the routine will also fix continuous
	variables, but the test still assumes a delta of 1.0.
    */
    virtual int reducedCostFix(double gap, bool justInteger=true) = 0 ;
  //@}
    
  /*! \name Methods to set variable type */
  //@{
    /** Set the index-th variable to be a continuous variable */
    virtual void setContinuous(int index) = 0;
    /** Set the index-th variable to be an integer variable */
    virtual void setInteger(int index) = 0;
    /** Set the variables listed in indices (which is of length len) to be
	continuous variables */
    virtual void setContinuous(const int* indices, int len) = 0 ;
    /** Set the variables listed in indices (which is of length len) to be
	integer variables */
    virtual void setInteger(const int* indices, int len) = 0 ;
  //@}

    /*! \brief Data type for name vectors. */
    typedef std::vector<std::string> OsiNameVec ;

  /*! \name Methods for row and column names

    Osi defines three name management disciplines: `auto names' (0), `lazy
    names' (1), and `full names' (2). See the description of
    #OsiNameDiscipline for details. Changing the name discipline (via
    setIntParam()) will not automatically add or remove name information,
    but setting the discipline to auto will make existing information
    inaccessible until the discipline is reset to lazy or full.

    By definition, a row index of getNumRows() (<i>i.e.</i>, one larger than
    the largest valid row index) refers to the objective function.

    OSI users and implementors: While the OSI base class can define an
    interface and provide rudimentary support, use of names really depends
    on support by the OsiXXX class to ensure that names are managed
    correctly.  If an OsiXXX class does not support names, it should return
    false for calls to getIntParam() or setIntParam() that reference
    #OsiNameDiscipline.
  */
  //@{

    /*! \brief Generate a standard name of the form Rnnnnnnn or Cnnnnnnn
    
      Set \p rc to 'r' for a row name, 'c' for a column name.
      The `nnnnnnn' part is generated from ndx and will contain 7 digits
      by default, padded with zeros if necessary. As a special case,
      ndx = getNumRows() is interpreted as a request for the name of the
      objective function. OBJECTIVE is returned, truncated to digits+1
      characters to match the row and column names.
    */
    virtual std::string dfltRowColName(char rc,
				 int ndx, unsigned digits = 7) const = 0 ;

    /*! \brief Return the name of the objective function */

  virtual std::string getObjName (unsigned maxLen = static_cast<unsigned>(std::string::npos)) const = 0 ;

    /*! \brief Set the name of the objective function */

    virtual void setObjName (std::string name) = 0 ;

    /*! \brief Return the name of the row.

      The routine will <i>always</i> return some name, regardless of the name
      discipline or the level of support by an OsiXXX derived class. Use
      maxLen to limit the length.
    */
    virtual std::string getRowName(int rowIndex,
				   unsigned maxLen = static_cast<unsigned>(std::string::npos)) const = 0 ;

    /*! \brief Return a pointer to a vector of row names

      If the name discipline (#OsiNameDiscipline) is auto, the return value
      will be a vector of length zero. If the name discipline is lazy, the
      vector will contain only names supplied by the client and will be no
      larger than needed to hold those names; entries not supplied will be
      null strings. In particular, the objective name is <i>not</i>
      included in the vector for lazy names. If the name discipline is
      full, the vector will have getNumRows() names, either supplied or
      generated, plus one additional entry for the objective name.
    */
    virtual const OsiNameVec &getRowNames() = 0 ;

    /*! \brief Set a row name

      Quietly does nothing if the name discipline (#OsiNameDiscipline) is
      auto. Quietly fails if the row index is invalid.
    */
    virtual void setRowName(int ndx, std::string name) = 0 ;

    /*! \brief Set multiple row names

      The run of len entries starting at srcNames[srcStart] are installed as
      row names starting at row index tgtStart. The base class implementation
      makes repeated calls to setRowName.
    */
    virtual void setRowNames(OsiNameVec &srcNames,
		     int srcStart, int len, int tgtStart) = 0 ;

    /*! \brief Delete len row names starting at index tgtStart

      The specified row names are removed and the remaining row names are
      copied down to close the gap.
    */
    virtual void deleteRowNames(int tgtStart, int len) = 0 ;
  
    /*! \brief Return the name of the column

      The routine will <i>always</i> return some name, regardless of the name
      discipline or the level of support by an OsiXXX derived class. Use
      maxLen to limit the length.
    */
    virtual std::string getColName(int colIndex,
				   unsigned maxLen = static_cast<unsigned>(std::string::npos)) const = 0 ;

    /*! \brief Return a pointer to a vector of column names

      If the name discipline (#OsiNameDiscipline) is auto, the return value
      will be a vector of length zero. If the name discipline is lazy, the
      vector will contain only names supplied by the client and will be no
      larger than needed to hold those names; entries not supplied will be
      null strings. If the name discipline is full, the vector will have
      getNumCols() names, either supplied or generated.
    */
    virtual const OsiNameVec &getColNames() = 0 ;

    /*! \brief Set a column name

      Quietly does nothing if the name discipline (#OsiNameDiscipline) is
      auto. Quietly fails if the column index is invalid.
    */
    virtual void setColName(int ndx, std::string name) = 0 ;

    /*! \brief Set multiple column names

      The run of len entries starting at srcNames[srcStart] are installed as
      column names starting at column index tgtStart. The base class
      implementation makes repeated calls to setColName.
    */
    virtual void setColNames(OsiNameVec &srcNames,
		     int srcStart, int len, int tgtStart) = 0 ;

    /*! \brief Delete len column names starting at index tgtStart

      The specified column names are removed and the remaining column names
      are copied down to close the gap.
    */
    virtual void deleteColNames(int tgtStart, int len) = 0 ;
  

    /*! \brief Set row and column names from a CoinMpsIO object.
    
      Also sets the name of the objective function. If the name discipline
      is auto, you get what you asked for. This routine does not use
      setRowName or setColName.
    */
    virtual void setRowColNames(const CoinMpsIO &mps) = 0 ;

    /*! \brief Set row and column names from a CoinModel object.

      If the name discipline is auto, you get what you asked for.
      This routine does not use setRowName or setColName.
    */
    virtual void setRowColNames(CoinModel &mod) = 0 ;

    /*! \brief Set row and column names from a CoinLpIO object.

      Also sets the name of the objective function. If the name discipline is
      auto, you get what you asked for. This routine does not use setRowName
      or setColName.
    */
    virtual void setRowColNames(CoinLpIO &mod) = 0 ;

  //@}

  /*! \name Methods to modify the constraint system.

     Note that new columns are added as continuous variables.
  */
  //@{

    /** Add a column (primal variable) to the problem. */
    virtual void addCol(const CoinPackedVectorBase& vec,
			const double collb, const double colub,   
			const double obj) = 0 ;

    /*! \brief Add a named column (primal variable) to the problem.

      The default implementation adds the column, then changes the name. This
      can surely be made more efficient within an OsiXXX class.
    */
    virtual void addCol(const CoinPackedVectorBase& vec,
			const double collb, const double colub,   
			const double obj, std::string name) = 0 ;

    /** Add a column (primal variable) to the problem. */
    virtual void addCol(int numberElements,
			const int* rows, const double* elements,
			const double collb, const double colub,   
			const double obj) = 0 ;

    /*! \brief Add a named column (primal variable) to the problem.

      The default implementation adds the column, then changes the name. This
      can surely be made more efficient within an OsiXXX class.
    */
    virtual void addCol(int numberElements,
			const int* rows, const double* elements,
			const double collb, const double colub,   
			const double obj, std::string name) = 0 ;

    /** Add a set of columns (primal variables) to the problem.
    
      The default implementation simply makes repeated calls to
      addCol().
    */
    virtual void addCols(const int numcols,
			 const CoinPackedVectorBase * const * cols,
			 const double* collb, const double* colub,   
			 const double* obj) = 0 ;

    /** Add a set of columns (primal variables) to the problem.
    
      The default implementation simply makes repeated calls to
      addCol().
    */
    virtual void addCols(const int numcols, const int* columnStarts,
			 const int* rows, const double* elements,
			 const double* collb, const double* colub,   
			 const double* obj) = 0 ;

    /// Add columns using a CoinBuild object
    virtual void addCols(const CoinBuild & buildObject) = 0 ;

    /** Add columns from a model object.  returns
       -1 if object in bad state (i.e. has row information)
       otherwise number of errors
       modelObject non const as can be regularized as part of build
    */
    virtual int addCols(CoinModel & modelObject) = 0 ;

    /** \brief Remove a set of columns (primal variables) from the
	       problem.

      The solver interface for a basis-oriented solver will maintain valid
      warm start information if all deleted variables are nonbasic.
    */
    virtual void deleteCols(const int num, const int * colIndices) = 0;
  
    /*! \brief Add a row (constraint) to the problem. */
    virtual void addRow(const CoinPackedVectorBase& vec,
			const double rowlb, const double rowub) = 0;

    /*! \brief Add a named row (constraint) to the problem.
    
      The default implementation adds the row, then changes the name. This
      can surely be made more efficient within an OsiXXX class.
    */
    virtual void addRow(const CoinPackedVectorBase& vec,
			const double rowlb, const double rowub,
			std::string name) = 0 ;

    /*! \brief Add a row (constraint) to the problem. */
    virtual void addRow(const CoinPackedVectorBase& vec,
			const char rowsen, const double rowrhs,   
			const double rowrng) = 0;

    /*! \brief Add a named row (constraint) to the problem.

      The default implementation adds the row, then changes the name. This
      can surely be made more efficient within an OsiXXX class.
    */
    virtual void addRow(const CoinPackedVectorBase& vec,
			const char rowsen, const double rowrhs,   
			const double rowrng, std::string name) = 0 ;

    /*! Add a row (constraint) to the problem.
    
      Converts to addRow(CoinPackedVectorBase&,const double,const double).
    */
    virtual void addRow(int numberElements,
			const int *columns, const double *element,
			const double rowlb, const double rowub) = 0 ;

    /*! Add a set of rows (constraints) to the problem.
    
      The default implementation simply makes repeated calls to
      addRow().
    */
    virtual void addRows(const int numrows,
			 const CoinPackedVectorBase * const * rows,
			 const double* rowlb, const double* rowub) = 0 ;

    /** Add a set of rows (constraints) to the problem.
    
      The default implementation simply makes repeated calls to
      addRow().
    */
    virtual void addRows(const int numrows,
			 const CoinPackedVectorBase * const * rows,
			 const char* rowsen, const double* rowrhs,   
			 const double* rowrng) = 0 ;

    /** Add a set of rows (constraints) to the problem.
    
      The default implementation simply makes repeated calls to
      addRow().
    */
    virtual void addRows(const int numrows, const int *rowStarts,
			 const int *columns, const double *element,
			 const double *rowlb, const double *rowub) = 0 ;

    /// Add rows using a CoinBuild object
    virtual void addRows(const CoinBuild &buildObject) = 0 ;

    /*! Add rows from a CoinModel object.
    
      Returns -1 if the object is in the wrong state (<i>i.e.</i>, has
      column-major information), otherwise the number of errors.

      The modelObject is not const as it can be regularized as part of
      the build.
    */
    virtual int addRows(CoinModel &modelObject) = 0 ;

    /** \brief Delete a set of rows (constraints) from the problem.

      The solver interface for a basis-oriented solver will maintain valid
      warm start information if all deleted rows are loose.
    */
    virtual void deleteRows(const int num, const int * rowIndices) = 0;

    /** \brief Replace the constraint matrix

      I (JJF) am getting annoyed because I can't just replace a matrix.
      The default behavior of this is do nothing so only use where that would
      not matter, e.g. strengthening a matrix for MIP.
    */
    virtual void replaceMatrixOptional(const CoinPackedMatrix & ) = 0 ;

    /** \brief Replace the constraint matrix
    
      And if it does matter (not used at present)
    */
    virtual void replaceMatrix(const CoinPackedMatrix & ) = 0 ;

    /** \brief Save a copy of the base model
    
      If solver wants it can save a copy of "base" (continuous) model here.
    */
    virtual void saveBaseModel() = 0 ;

    /** \brief Reduce the constraint system to the specified number of
    	       constraints.

       If solver wants it can restore a copy of "base" (continuous) model
       here.

       \note
       The name is somewhat misleading. Implementors should consider
       the opportunity to optimise behaviour in the common case where
       \p numberRows is exactly the number of original constraints. Do not,
       however, neglect the possibility that \p numberRows does not equal
       the number of original constraints.
     */
    virtual void restoreBaseModel(int numberRows) = 0 ;
    //-----------------------------------------------------------------------
    /** Apply a collection of cuts.

	Only cuts which have an <code>effectiveness >= effectivenessLb</code>
	are applied.
	<ul>
	  <li> ReturnCode.getNumineffective() -- number of cuts which were
	       not applied because they had an
	       <code>effectiveness < effectivenessLb</code>
	  <li> ReturnCode.getNuminconsistent() -- number of invalid cuts
	  <li> ReturnCode.getNuminconsistentWrtIntegerModel() -- number of
	       cuts that are invalid with respect to this integer model
	  <li> ReturnCode.getNuminfeasible() -- number of cuts that would
	       make this integer model infeasible
	  <li> ReturnCode.getNumApplied() -- number of integer cuts which
	       were applied to the integer model
	  <li> cs.size() == getNumineffective() +
			    getNuminconsistent() +
			    getNuminconsistentWrtIntegerModel() +
			    getNuminfeasible() +
			    getNumApplied()
	</ul>
    */
    virtual ApplyCutsReturnCode applyCuts(const OsiCuts & cs,
					  double effectivenessLb = 0.0)
    { 
      return (applyCutsPrivate(cs,effectivenessLb)) ;
    }

    /** Apply a collection of row cuts which are all effective.
	applyCuts seems to do one at a time which seems inefficient.
	Would be even more efficient to pass an array of pointers.
    */
    virtual void applyRowCuts(int numberCuts, const OsiRowCut * cuts) = 0 ;

    /** Apply a collection of row cuts which are all effective.
	This is passed in as an array of pointers.
    */
    virtual void applyRowCuts(int numberCuts, const OsiRowCut ** cuts) = 0 ;

    /// Deletes branching information before columns deleted
    virtual void deleteBranchingInfo(int numberDeleted, const int * which) = 0 ;

  //@}

  /*! \name Methods for problem input and output */
  //@{
    /*! \brief Load in a problem by copying the arguments. The constraints on
	    the rows are given by lower and upper bounds.
	
	If a pointer is 0 then the following values are the default:
        <ul>
          <li> <code>colub</code>: all columns have upper bound infinity
          <li> <code>collb</code>: all columns have lower bound 0 
          <li> <code>rowub</code>: all rows have upper bound infinity
          <li> <code>rowlb</code>: all rows have lower bound -infinity
	  <li> <code>obj</code>: all variables have 0 objective coefficient
        </ul>

	Note that the default values for rowub and rowlb produce the
	constraint -infty <= ax <= infty. This is probably not what you want.
    */
    virtual void loadProblem (const CoinPackedMatrix& matrix,
			      const double* collb, const double* colub,   
			      const double* obj,
			      const double* rowlb, const double* rowub) = 0;
			    
    /*! \brief Load in a problem by assuming ownership of the arguments.
	    The constraints on the rows are given by lower and upper bounds.

	For default argument values see the matching loadProblem method.

	\warning
	The arguments passed to this method will be freed using the
	C++ <code>delete</code> and <code>delete[]</code> functions. 
    */
    virtual void assignProblem (CoinPackedMatrix*& matrix,
			        double*& collb, double*& colub, double*& obj,
			        double*& rowlb, double*& rowub) = 0;

    /*! \brief Load in a problem by copying the arguments.
	    The constraints on the rows are given by sense/rhs/range triplets.
	    
	If a pointer is 0 then the following values are the default:
	<ul>
          <li> <code>colub</code>: all columns have upper bound infinity
          <li> <code>collb</code>: all columns have lower bound 0 
	  <li> <code>obj</code>: all variables have 0 objective coefficient
          <li> <code>rowsen</code>: all rows are >=
          <li> <code>rowrhs</code>: all right hand sides are 0
          <li> <code>rowrng</code>: 0 for the ranged rows
        </ul>

	Note that the default values for rowsen, rowrhs, and rowrng produce the
	constraint ax >= 0.
    */
    virtual void loadProblem (const CoinPackedMatrix& matrix,
			      const double* collb, const double* colub,
			      const double* obj,
			      const char* rowsen, const double* rowrhs,   
			      const double* rowrng) = 0;

    /*! \brief Load in a problem by assuming ownership of the arguments.
	    The constraints on the rows are given by sense/rhs/range triplets.
	
	For default argument values see the matching loadProblem method.

	\warning
	The arguments passed to this method will be freed using the
	C++ <code>delete</code> and <code>delete[]</code> functions. 
    */
    virtual void assignProblem (CoinPackedMatrix*& matrix,
			        double*& collb, double*& colub, double*& obj,
			        char*& rowsen, double*& rowrhs,
			        double*& rowrng) = 0;

    /*! \brief Load in a problem by copying the arguments. The constraint
	    matrix is is specified with standard column-major
	    column starts / row indices / coefficients vectors. 
	    The constraints on the rows are given by lower and upper bounds.
    
      The matrix vectors must be gap-free. Note that <code>start</code> must
      have <code>numcols+1</code> entries so that the length of the last column
      can be calculated as <code>start[numcols]-start[numcols-1]</code>.

      See the previous loadProblem method using rowlb and rowub for default
      argument values.
    */
    virtual void loadProblem (const int numcols, const int numrows,
			      const CoinBigIndex * start, const int* index,
			      const double* value,
			      const double* collb, const double* colub,   
			      const double* obj,
			      const double* rowlb, const double* rowub) = 0;

    /*! \brief Load in a problem by copying the arguments. The constraint
	    matrix is is specified with standard column-major
	    column starts / row indices / coefficients vectors. 
	    The constraints on the rows are given by sense/rhs/range triplets.
    
      The matrix vectors must be gap-free. Note that <code>start</code> must
      have <code>numcols+1</code> entries so that the length of the last column
      can be calculated as <code>start[numcols]-start[numcols-1]</code>.

      See the previous loadProblem method using sense/rhs/range for default
      argument values.
    */
    virtual void loadProblem (const int numcols, const int numrows,
			      const CoinBigIndex * start, const int* index,
			      const double* value,
			      const double* collb, const double* colub,   
			      const double* obj,
			      const char* rowsen, const double* rowrhs,   
  			      const double* rowrng) = 0;

    /*! \brief Load a model from a CoinModel object. Return the number of
	    errors encountered.

      The modelObject parameter cannot be const as it may be changed as part
      of process. If keepSolution is true will try and keep warmStart.
    */
    virtual int loadFromCoinModel (CoinModel & modelObject,
				   bool keepSolution=false) = 0 ;

    /*! \brief Read a problem in MPS format from the given filename.
    
      The default implementation uses CoinMpsIO::readMps() to read
      the MPS file and returns the number of errors encountered.
    */
    virtual int readMps (const char *filename,
			 const char *extension = "mps") = 0 ;

    /*! \brief Read a problem in MPS format from the given full filename.
    
      This uses CoinMpsIO::readMps() to read the MPS file and returns the
      number of errors encountered. It also may return an array of set
      information
    */
    virtual int readMps (const char *filename, const char*extension,
			int & numberSets, CoinSet ** & sets) = 0 ;

    /*! \brief Read a problem in GMPL format from the given filenames.
    
      The default implementation uses CoinMpsIO::readGMPL(). This capability
      is available only if the third-party package Glpk is installed.
    */
    virtual int readGMPL (const char *filename, const char *dataname=nullptr) = 0 ;

    /*! \brief Write the problem in MPS format to the specified file.

      If objSense is non-zero, a value of -1.0 causes the problem to be
      written with a maximization objective; +1.0 forces a minimization
      objective. If objSense is zero, the choice is left to the implementation.
    */
    virtual void writeMps (const char *filename,
			   const char *extension = "mps",
			   double objSense=0.0) const = 0;

    /*! \brief Write the problem in MPS format to the specified file with
	    more control over the output.

	Row and column names may be null.
	formatType is
	<ul>
	  <li> 0 - normal
	  <li> 1 - extra accuracy 
	  <li> 2 - IEEE hex
	</ul>

	Returns non-zero on I/O error
    */
    virtual int writeMpsNative (const char *filename, 
		        const char ** rowNames, const char ** columnNames,
		        int formatType=0,int numberAcross=2,
		        double objSense=0.0, int numberSOS=0,
		        const CoinSet * setInfo=nullptr) const = 0 ;

/***********************************************************************/
// Lp files 

  /** Write the problem into an Lp file of the given filename with the 
      specified extension.
      Coefficients with value less than epsilon away from an integer value
      are written as integers.
      Write at most numberAcross monomials on a line.
      Write non integer numbers with decimals digits after the decimal point.

      The written problem is always a minimization problem.
      If the current problem is a maximization problem, the 
      intended objective function for the written problem is the current
      objective function multiplied by -1. If the current problem is a
      minimization problem, the intended objective function for the
      written problem is the current objective function.
      If objSense < 0, the intended objective function is multiplied by -1
      before writing the problem. It is left unchanged otherwise.

      Write objective function name and constraint names if useRowNames is 
      true. This version calls writeLpNative().
  */
  virtual void writeLp(const char *filename,
               const char *extension = "lp",
                double epsilon = 1e-5,
                int numberAcross = 10,
                int decimals = 5,
                double objSense = 0.0,
	        bool useRowNames = true) const = 0 ;

  /** Write the problem into the file pointed to by the parameter fp. 
      Other parameters are similar to 
      those of writeLp() with first parameter filename.
  */
  virtual void writeLp(FILE *fp,
                double epsilon = 1e-5,
                int numberAcross = 10,
                int decimals = 5,
                double objSense = 0.0,
	        bool useRowNames = true) const = 0 ;

  /** Write the problem into an Lp file. Parameters are similar to 
      those of writeLp(), but in addition row names and column names
      may be given. 

      Parameter rowNames may be null, in which case default row names 
      are used. If rowNames is not null, it must have exactly one entry
      per row in the problem and one additional
      entry (rowNames[getNumRows()] with the objective function name.
      These getNumRows()+1 entries must be distinct. If this is not the 
      case, default row names
      are used. In addition, format restrictions are imposed on names
      (see CoinLpIO::is_invalid_name() for details).

      Similar remarks can be made for the parameter columnNames which
      must either be null or have exactly getNumCols() distinct entries.

      Write objective function name and constraint names if 
      useRowNames is true. */
  virtual int writeLpNative(const char *filename,
		    char const * const * const rowNames,
		    char const * const * const columnNames,
		    const double epsilon = 1.0e-5,
                    const int numberAcross = 10,
                    const int decimals = 5,
                    const double objSense = 0.0,
		    const bool useRowNames = true) const = 0 ;

  /** Write the problem into the file pointed to by the parameter fp. 
      Other parameters are similar to 
      those of writeLpNative() with first parameter filename.
  */
  virtual int writeLpNative(FILE *fp,
		    char const * const * const rowNames,
		    char const * const * const columnNames,
		    const double epsilon = 1.0e-5,
                    const int numberAcross = 10,
                    const int decimals = 5,
                    const double objSense = 0.0,
		    const bool useRowNames = true) const = 0 ;

  /// Read file in LP format from file with name filename. 
  /// See class CoinLpIO for description of this format.
  virtual int readLp(const char *filename, const double epsilon = 1e-5) = 0 ;

  /// Read file in LP format from the file pointed to by fp. 
  /// See class CoinLpIO for description of this format.
  virtual int readLp(FILE *fp, const double epsilon = 1e-5) = 0 ;

  //@}

  /*! \name Setting/Accessing application data */
  //@{
    /** Set application data.

	This is a pointer that the application can store into and
	retrieve from the solver interface.
	This field is available for the application to optionally
	define and use.
    */
    virtual void setApplicationData (void * appData) = 0 ;
    /** Create a clone of an Auxiliary Information object.
        The base class just stores an application data pointer
        but can be more general.  Application data pointer is
        designed for one user while this can be extended to cope
        with more general extensions.
    */
    virtual void setAuxiliaryInfo(OsiAuxInfo * auxiliaryInfo) = 0 ;

    /// Get application data
    virtual void * getApplicationData() const = 0 ;
    /// Get pointer to auxiliary info object
    virtual OsiAuxInfo * getAuxiliaryInfo() const = 0 ;
  //@}

  /*! \name Message handling
  
    See the COIN library documentation for additional information about
    COIN message facilities.
  
  */
  //@{
  /** Pass in a message handler
  
    It is the client's responsibility to destroy a message handler installed
    by this routine; it will not be destroyed when the solver interface is
    destroyed. 
  */
  virtual void passInMessageHandler(CoinMessageHandler * handler) = 0 ;
  /// Set language
  virtual void setLanguage(CoinMessages::Language language) = 0 ;
  /// Return a pointer to the current message handler
  virtual CoinMessageHandler * messageHandler() const = 0 ;
  /// Return the current set of messages
  virtual CoinMessages messages() = 0 ;
  /// Return a pointer to the current set of messages
  virtual CoinMessages * messagesPointer()  = 0 ;
  /// Return true if default handler
  virtual bool defaultHandler() const = 0 ;
  //@}

  /*! \name Methods for dealing with discontinuities other than integers.
  
     Osi should be able to know about SOS and other types.  This is an optional
     section where such information can be stored.
  */
  //@{
    /** \brief Identify integer variables and create corresponding objects.
  
      Record integer variables and create an OsiSimpleInteger object for each
      one.  All existing OsiSimpleInteger objects will be destroyed.
      If justCount then no objects created and we just store numberIntegers_
    */
    virtual void findIntegers(bool justCount) = 0 ;

    /** \brief Identify integer variables and SOS and create corresponding
	       objects.
  
      Record integer variables and create an OsiSimpleInteger object for each
      one.  All existing OsiSimpleInteger objects will be destroyed.
      If the solver supports SOS then do the same for SOS.

      If justCount then no objects created and we just store numberIntegers_
      Returns number of SOS
    */
    virtual int findIntegersAndSOS(bool justCount) = 0 ;

    /// Get the number of objects
    virtual int numberObjects() const = 0 ;

    /// Set the number of objects
    virtual void setNumberObjects(int number) = 0 ;

    /// Get the array of objects
    virtual OsiObject ** objects() const = 0 ;

    /// Get the specified object
    const virtual OsiObject * object(int which) const = 0 ;

    /// Get the specified object
    virtual OsiObject * modifiableObject(int which) const = 0 ;

    /// Delete all object information
    virtual void deleteObjects() = 0 ;

    /** Add in object information.
  
      Objects are cloned; the owner can delete the originals.
    */
    virtual void addObjects(int numberObjects, OsiObject ** objects) = 0 ;
    /*! \brief Use current solution to set bounds so current integer feasible
	       solution will stay feasible.

        Only feasible bounds will be used, even if the current solution is
	outside bounds.  The amount of such violation will be returned (and
	if small can be ignored)
    */
    virtual double forceFeasible() = 0 ;
  //@}

  /*! \name Methods related to testing generated cuts
  
    See the documentation for OsiRowCutDebugger for additional details.
  */
  //@{
    /*! \brief Activate the row cut debugger.

      If \p modelName is in the set of known models then all cuts are
      checked to see that they do NOT cut off the optimal solution known
      to the debugger.
    */
    virtual void activateRowCutDebugger (const char *modelName) = 0 ;

    /*! \brief Activate the row cut debugger using a full solution array.


      Activate the debugger for a model not included in the debugger's
      internal database.  Cuts will be checked to see that they do NOT
      cut off the given solution.

      \p solution must be a full solution vector, but only the integer
      variables need to be correct. The debugger will fill in the continuous
      variables by solving an lp relaxation with the integer variables
      fixed as specified. If the given values for the continuous variables
      should be preserved, set \p keepContinuous to true.
    */
    virtual void activateRowCutDebugger(const double *solution,
    					bool enforceOptimality = true) = 0 ;

    /*! \brief Get the row cut debugger provided the solution known to the
    	       debugger is within the feasible region held in the solver.

      If there is a row cut debugger object associated with model AND if
      the solution known to the debugger is within the solver's current
      feasible region (i.e., the column bounds held in the solver are
      compatible with the known solution) then a pointer to the debugger
      is returned which may be used to test validity of cuts.

      Otherwise nullptr is returned
    */
    virtual const OsiRowCutDebugger *getRowCutDebugger() const = 0 ;

    /*! \brief Get the row cut debugger object

      Return the row cut debugger object if it exists. One common usage of
      this method is to obtain a debugger object in order to execute
      OsiRowCutDebugger::redoSolution (so that the stored solution is again
      compatible with the problem held in the solver).
    */
    virtual OsiRowCutDebugger * getRowCutDebuggerAlways() const = 0 ;
  //@} 

  /*! \name OSI simplex interface

    Methods for an advanced interface to a simplex solver. The interface
    comprises two groups of methods. Group 1 contains methods for tableau
    access. Group 2 contains methods for dictating individual simplex pivots.
  */
  //@{

  /*! \brief Return the simplex implementation level.
  
      The return codes are:
      - 0: the simplex interface is not implemented.
      - 1: the Group 1 (tableau access) methods are implemented.
      - 2: the Group 2 (pivoting) methods are implemented

      The codes are cumulative - a solver which implements Group 2 also
      implements Group 1.
  */
  virtual int canDoSimplexInterface() const = 0 ;
  //@}

  /*! \name OSI simplex interface Group 1
      \brief Tableau access methods.

      This group of methods provides access to rows and columns of the basis
      inverse and to rows and columns of the tableau.
  */
  //@{

  /*! \brief Prepare the solver for the use of tableau access methods.
  
    Prepares the solver for the use of the tableau access methods, if
    any such preparation is required.

    The \c const attribute is required due to the places this method
    may be called (e.g., within CglCutGenerator::generateCuts()).
  */
  virtual void enableFactorization() const = 0 ;

  /*! \brief Undo the effects of #enableFactorization. */
  virtual void disableFactorization() const = 0 ;

  /*! \brief Check if an optimal basis is available.
  
    Returns true if the problem has been solved to optimality and a
    basis is available. This should be used to see if the tableau access
    operations are possible and meaningful.
    
    \note
    Implementors please note that this method may be called
    before #enableFactorization.
  */
  virtual bool basisIsAvailable() const = 0 ;

  /// Synonym for #basisIsAvailable
  inline bool optimalBasisIsAvailable() const { return basisIsAvailable() ; }

  /*! \brief Retrieve status information for column and row variables.
 
    This method returns status as integer codes:
    <ul>
      <li> 0: free
      <li> 1: basic
      <li> 2: nonbasic at upper bound
      <li> 3: nonbasic at lower bound
    </ul>

    The #getWarmStart method provides essentially the same functionality
    for a simplex-oriented solver, but the implementation details are very
    different.

    \note
    Logical variables associated with rows are all assumed to have +1
    coefficients, so for a <= constraint the logical will be at lower
    bound if the constraint is tight.

    \note
    Implementors may choose to implement this method as a wrapper which
    converts a CoinWarmStartBasis to the requested representation.
  */
  virtual void getBasisStatus(int* cstat, int* rstat) const = 0 ;

  /*! \brief Set the status of column and row variables and update
             the basis factorization and solution.

    Status information should be coded as documented for #getBasisStatus.
    Returns 0 if all goes well, 1 if something goes wrong.

    This method differs from #setWarmStart in the format of the input
    and in its immediate effect. Think of it as #setWarmStart immediately
    followed by #resolve, but no pivots are allowed.

    \note
    Implementors may choose to implement this method as a wrapper that calls
    #setWarmStart and #resolve if the no pivot requirement can be satisfied.
  */
  virtual int setBasisStatus(const int* cstat, const int* rstat) = 0 ;

  /*! \brief Calculate duals and reduced costs for the given objective
	     coefficients.

    The solver's objective coefficient vector is not changed.
  */
  virtual void getReducedGradient(double* columnReducedCosts, 
				  double* duals, const double* c) const = 0 ;

  /*! \brief Get a row of the tableau
  
    If \p slack is not null, it will be loaded with the coefficients for
    the artificial (logical) variables (i.e., the row of the basis inverse).
  */
  virtual void getBInvARow(int row, double* z, double* slack = nullptr) const = 0 ;

  /*! \brief Get a row of the basis inverse */
  virtual void getBInvRow(int row, double* z) const = 0 ;

  /*! \brief Get a column of the tableau */
  virtual void getBInvACol(int col, double* vec) const = 0 ;

  /*! \brief Get a column of the basis inverse */
  virtual void getBInvCol(int col, double* vec) const = 0 ;

  /*! \brief Get indices of basic variables
  
    If the logical (artificial) for row i is basic, the index should be coded
    as (#getNumCols + i).
    The order of indices must match the order of elements in the vectors
    returned by #getBInvACol and #getBInvCol.
  */
  virtual void getBasics(int* index) const = 0 ;

  //@}

  /*! \name OSI simplex interface Group 2
      \brief Pivoting methods

      This group of methods provides for control of individual pivots by a
      simplex solver.
  */
  //@{

  /**Enables normal operation of subsequent functions.
     This method is supposed to ensure that all typical things (like
     reduced costs, etc.) are updated when individual pivots are executed
     and can be queried by other methods.  says whether will be
     doing primal or dual
  */
  virtual void enableSimplexInterface(bool doingPrimal) = 0 ;

  ///Undo whatever setting changes the above method had to make
  virtual void disableSimplexInterface() = 0 ;

  /** Perform a pivot by substituting a colIn for colOut in the basis. 
     The status of the leaving variable is given in outStatus. Where
     1 is to upper bound, -1 to lower bound
     Return code was undefined - now for OsiClp is 0 for okay,
     1 if inaccuracy forced re-factorization (should be okay) and
     -1 for singular factorization
  */
  virtual int pivot(int colIn, int colOut, int outStatus) = 0 ;

  /** Obtain a result of the primal pivot 
      Outputs: colOut -- leaving column, outStatus -- its status,
      t -- step size, and, if dx!=nullptr, *dx -- primal ray direction.
      Inputs: colIn -- entering column, sign -- direction of its change (+/-1).
      Both for colIn and colOut, artificial variables are index by
      the negative of the row index minus 1.
      Return code (for now): 0 -- leaving variable found, 
      -1 -- everything else?
      Clearly, more informative set of return values is required 
      Primal and dual solutions are updated
  */
  virtual int primalPivotResult(int colIn, int sign, 
				int& colOut, int& outStatus, 
				double& t, CoinPackedVector* dx) = 0 ;

  /** Obtain a result of the dual pivot (similar to the previous method)
      Differences: entering variable and a sign of its change are now
      the outputs, the leaving variable and its statuts -- the inputs
      If dx!=nullptr, then *dx contains dual ray
      Return code: same
  */
  virtual int dualPivotResult(int& colIn, int& sign, 
			      int colOut, int outStatus, 
			      double& t, CoinPackedVector* dx) = 0 ;
  //@}

  protected:

  /*! \brief We need this for covariant return */
  virtual ApplyCutsReturnCode applyCutsPrivate(const OsiCuts & cs,
					  double effectivenessLb = 0.0) = 0 ;

  private:

  /// Intermediate variable for returning the API name
  const char *apiName_ ;
};

}  // end namespace Osi2

#endif
