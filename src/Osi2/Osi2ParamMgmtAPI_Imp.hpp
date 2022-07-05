/*
  Copyright 2019 Lou Hafer, Matt Saltzman
  This code is licensed under the terms of the Eclipse Public License (EPL)

  $Id$
*/
/*! \file Osi2ParamMgmtAPI_Imp.hpp

  Defines an implementation of Osi2::ParamMgmtAPI, an API to provide some
  uniformity and syntactic sugar for parameter management.
*/


#ifndef Osi2ParamMgmtAPI_Imp_HPP
#define Osi2ParamMgmtAPI_Imp_HPP

#include "Osi2API.hpp"
#include "Osi2PluginManager.hpp"

#include "Osi2ParamMgmtAPI.hpp"
#include "Osi2ParamMgmtAPIMessages.hpp"

namespace Osi2 {

/*! \brief Osi2 ParamMgmtAPI implementation

  This class implements the Osi2::ParamMgmt API.
*/
class OSI2LIB_EXPORT ParamMgmtAPI_Imp : public ParamMgmtAPI {

public:

  /// ParamMgmtAPI_Imp implements ParamMgmtAPI
  inline static const char *getAPIIDString ()
  { return (ParamMgmtAPI::getAPIIDString()) ; }

/*! \name Constructors and Destructors */
//@{
  /// Default constructor
  ParamMgmtAPI_Imp() ;
  /// Registration constructor
  ParamMgmtAPI_Imp(std::string name) ;
  /// Default copy constructor
  ParamMgmtAPI_Imp(const ParamMgmtAPI_Imp &original) ;
  /// Assignment
  ParamMgmtAPI_Imp &operator=(const ParamMgmtAPI_Imp &rhs) ;

  /// Create (factory constructor)
  ParamMgmtAPI *create() ;
  /// Clone (factory copy constructor)
  ParamMgmtAPI *clone() ;

  /// Destructor
  ~ParamMgmtAPI_Imp() ;
//@}

/*! \brief Parameter handling methods.

  Methods to enroll / remove the parameter set for an object, and to set
  and get the values of parameters.
*/
//@{

  /// Enroll an object and its parameters with the manager
  bool enroll(std::string ident, Osi2::API *object) ;

  /// Remove an object and its parameters from the manager
  bool remove(std::string ident) ;

  /// Set a parameter
  bool set(std::string ident, std::string param, const void *blob) ;

  /// Get a parameter
  bool get(std::string ident, std::string param, void *blob) ;

//@}


/*! \name ParamMgmt API control methods

  Miscellaneous methods that control the behaviour of a ParamMgmtAPI object.
*/
//@{

  /// Set the log (verbosity) level
  inline void setLogLvl(int logLvl) {
    logLvl_ = logLvl ;
    msgHandler_->setLogLevel(logLvl_) ;
  }

  /// Get the log (verbosity) level
  inline int getLogLvl() const { return (logLvl_) ; }

  /*! \brief Set the message handler

    Replaces the current message handler. If newHandler is null, the existing
    handler is replaced with a default handler. It is the responsibility of
    the client to destroy any handler it supplies. The plugin manager takes
    responsibility for a default handler.
  */
  void setMsgHandler(CoinMessageHandler *newHandler) ;

  /// Get the message handler
  inline CoinMessageHandler *getMsgHandler() const {
      return (msgHandler_) ;
  }

  /// Report owner of message handler (false if owned by client).
  inline bool dfltHandler () const {
      return (dfltHandler_) ;
  }

//@}

private:

  /*! \brief Parameter management structures */
  //@{
  
  /// Entry in the #index used to manage enrolled objects
  struct ObjData {
    /// Enrolled object
    Osi2::API *enrolledObject_ ;
    /// Parameter back-end object to manipulate parameters of #enrolledObject
    ParamBEAPI *paramHandler_ ;
    /// Parameter names exported by #enrolledObject_
    std::vector<const char *> paramNames_ ;
  } ;

  /// Typedef for the map used to manage enrolled objects
  typedef std::map<std::string,ObjData> IndexMap ;

  /// Index map to track enrolled objects
  IndexMap indexMap_ ;

  /// Indicator; false if the message handler belongs to the client
  bool dfltHandler_ ;
  /// Message handler
  CoinMessageHandler *msgHandler_ ;
  /// Messages
  CoinMessages msgs_ ;
  /// Log (verbosity) level
  int logLvl_ ;

} ;

} // namespace Osi2 ;

#endif
