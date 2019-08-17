#ifndef Osi2RunParamsAPI_HPP
# define Osi2RunParamsAPI_HPP

/*! \file Osi2RunParamsAPI.hpp

  Provides an %API that supports a single run-time parameter object which can
  hold a collection of parameters to be applied to an object.
*/

#include <string>
#include <vector>

#include "Osi2API.hpp"

namespace Osi2 {

/*! \brief Provide a `run-time' parameter object

  This %API provides an object that can hold `run-time 'parameters for
  an object.

  The original idea was to provide a way to set up alternate sets of parameter
  values that could be applied together to a solver prior to invoking solve,
  to make it easy to set up and test alternate parameter settings.

  The concept, however, is fairly generic: Provide an object that can
  be loaded with a set of parameters that are to be applied together.
  RunParamsAPI provides a set of methods that can be used to load parameters
  into a RunParamsAPI object, and a set of methods that can be used to
  return the set of parameters. Parameters are grouped by data type.
  Supported data types are int, double, std::string, and void*
*/
class RunParamsAPI : public API {

public:

  /// Return the "ident" string for the RunParams %API
  inline static const char *getAPIIDString () { return ("RunParams") ; }

  /*! \name Methods to work with individual parameters */
  //@{

  /// Add an integer parameter, default value 0
  virtual void addIntParam(std::string id, int val = 0) = 0 ;
  /// Set the value of an existing integer parameter
  virtual void setIntParam(std::string id, int val) = 0 ;
  /// Retrieve the value of an integer parameter
  virtual int getIntParam(std::string id) const = 0 ;

  /// Add a double parameter, default value 0.0
  virtual void addDblParam(std::string id, double val = 0.0) = 0 ;
  /// Set the value of an existing double parameter
  virtual void setDblParam(std::string id, double val) = 0 ;
  /// Retrieve the value of a double parameter
  virtual double getDblParam(std::string id) const = 0 ;

  /// Add a string parameter, default value string()
  virtual void addStrParam(std::string id,
  			   std::string val = std::string()) = 0 ;
  /// Set the value of an existing string parameter
  virtual void setStrParam(std::string id, std::string val) = 0 ;
  /// Retrieve the value of a double parameter
  virtual std::string getStrParam(std::string id) const = 0 ;

  /// Add a void* parameter, default value nullptr
  virtual void addVoidParam(std::string id, void *val = nullptr) = 0 ;
  /// Set the value of an existing void* parameter
  virtual void setVoidParam(std::string id, void *val) = 0 ;
  /// Retrieve the value of a double parameter
  virtual void *getVoidParam(std::string id) const = 0 ;

  //@}

  /*! \name Methods to retrieve a parameter set */
  //@{

  /// Return the names of integer parameters
  virtual std::vector<std::string> getIntParamIds() const = 0 ;
  /// Return the names of double parameters
  virtual std::vector<std::string> getDblParamIds() const = 0 ;
  /// Return the names of string parameters
  virtual std::vector<std::string> getStrParamIds() const = 0 ;
  /// Return the names of void* parameters
  virtual std::vector<std::string> getVoidParamIds() const = 0 ;

  //@}

protected:

  /// Virtual destructor

  ~RunParamsAPI () { } ;

} ;

} // end namespace Osi2

#endif
