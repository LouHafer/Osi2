#ifndef Osi2RunParamsAPI_Imp_HPP
# define Osi2RunParamsAPI_Imp_HPP

/*! \file Osi2RunParamsAPI_Imp.hpp

  Provides a class to implement Osi2::RunParamsAPI.
*/

#include <string>
#include <map>

#include "Osi2RunParamsAPI.hpp"

namespace Osi2 {

/*! \brief Implementation class for RunParamsAPI

  Implements the abstract interface using wrappers around templated function
  invocation, to hide the templating from the user.
*/
class RunParamsAPI_Imp : public RunParamsAPI {

public:

  /// Registration constructor
  RunParamsAPI_Imp (std::string name) ;
  /// Default constructor
  RunParamsAPI_Imp () ;
  /// Destructor
  ~RunParamsAPI_Imp() ;

  template <class ParamType>
  using ParamMap = std::map<std::string,ParamType> ;

  /*! \name Methods to work with individual parameters */
  //@{

  void addIntParam(std::string id, int val = 0)
  { addParam<int>(intParams_,id,val) ; }
  void setIntParam(std::string id, int val)
  { setParam<int>(intParams_,id,val) ; }
  int getIntParam(std::string id) const
  { return (getParam<int>(intParams_,id)) ; }

  void addDblParam(std::string id, double val = 0.0)
  { addParam<double>(dblParams_,id,val) ; }
  void setDblParam(std::string id, double val)
  { setParam<double>(dblParams_,id,val) ; }
  double getDblParam(std::string id) const
  { return (getParam<double>(dblParams_,id)) ; }

  void addStrParam(std::string id, std::string val = std::string())
  { addParam<std::string>(strParams_,id,val) ; }
  void setStrParam(std::string id, std::string val)
  { setParam<std::string>(strParams_,id,val) ; }
  std::string getStrParam(std::string id) const
  { return (getParam<std::string>(strParams_,id)) ; }

  void addVoidParam(std::string id, void *val = nullptr)
  { addParam<void *>(voidParams_,id,val) ; }
  void setVoidParam(std::string id, void *val)
  { setParam<void *>(voidParams_,id,val) ; }
  void *getVoidParam(std::string id) const
  { return(getParam<void *>(voidParams_,id)) ; }

  //@}

  /*! \name Methods to retrieve parameter names */
  //@{

  // Return the names of integer parameters
  std::vector<std::string> getIntParamIds() const
  { return (getParamIds<int>(intParams_)) ; }
  // Return the names of double parameters
  std::vector<std::string> getDblParamIds() const
  { return (getParamIds<double>(dblParams_)) ; }
  // Return the names of std::string parameters
  std::vector<std::string> getStrParamIds() const
  { return (getParamIds<std::string>(strParams_)) ; }
  // Return the names of void* parameters
  std::vector<std::string> getVoidParamIds() const
  { return (getParamIds<void *>(voidParams_)) ; }

  //@}

private:

  /*! \name Constructors, destructor, etc. */
  //@{
  /// Copy constructor
  RunParamsAPI_Imp(RunParamsAPI_Imp &rhs) ;

  /// Virtual constructor
  RunParamsAPI *create() { return (new RunParamsAPI_Imp()) ; }
  /// Virtual copy constructor
  RunParamsAPI *clone(RunParamsAPI *rhs)
  { RunParamsAPI_Imp *tmp = dynamic_cast<RunParamsAPI_Imp *>(rhs) ;
    if (tmp != nullptr) {
      return (new RunParamsAPI_Imp(*tmp)) ;
    } else {
      return (new RunParamsAPI_Imp()) ;
    }
  }

  //@}

  /*! \brief Underlying template functions */
  //@{
  template <class ParamType>
  void addParam(ParamMap<ParamType> &paramMap, std::string id, ParamType val) ;
  template <class ParamType>
  void setParam(ParamMap<ParamType> &paramMap, std::string id, ParamType val) ;
  template <class ParamType>
  ParamType getParam(const ParamMap<ParamType> &paramMap,
  		     std::string id) const ;
  template <class ParamType>
  std::vector<std::string> getParamIds(const ParamMap<ParamType> &paramMap)
  				       		const ;
  //@}

  /*! \brief Parameter maps */
  //@{
  ParamMap<int> intParams_ ;
  ParamMap<double> dblParams_ ;
  ParamMap<std::string> strParams_ ;
  ParamMap<void *> voidParams_ ;
  //@}

} ;

/*
  Bring in the templated function definitions.
*/
#define Provide_RunParamsAPI_TemplateMethodDefs
#include "Osi2RunParamsAPI_Imp.cpp"
#undef Provide_RunParamsAPI_TemplateMethodDefs

}  // end namespace Osi2

#endif

