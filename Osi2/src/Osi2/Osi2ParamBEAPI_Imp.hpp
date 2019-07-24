/*
  Copyright 2019 Lou Hafer, Matt Saltzman, Kevin Poton
  This code is licensed under the terms of the Eclipse Public License (EPL)

  $Id$
*/
/*! \file Osi2ParamBEAPI_Imp.hpp

  This file defines a class ParamBEAPI_Imp that implements the ParamBEAPI
  interface.
*/


#ifndef Osi2ParamBEAPI_Imp_HPP
#define Osi2ParamBEAPI_Imp_HPP

#include <map>

#include "Osi2API.hpp"

#include "Osi2ParamBEAPI.hpp"
#include "Osi2ParamBEAPIMessages.hpp"

namespace Osi2 {

/*! \brief An implementation of the parameter management back end.

  This class implements the ParamBEAPI methods that must be supplied by
  an %API implementation to support parameter management. It is intended
  to be used by including a ParamBEAPI_Imp object as a member in the %API
  object that is to be managed.

  In particular, it can provide support for parameter management in an
  existing %API implementation with minimal modifications to the existing
  class. Given an implementation class SomeAPI_Imp that implements SomeAPI,
  first add a member to SomeAPI_Imp:
  \code{.cpp}
    ParamBEAPI_Imp<SomeAPI> paramHdlr ;
  \endcode

  Next, load the APIs supported by SomeAPI_Imp into the handler, along with a
  pointer to the object that supports the %API:
  \code{.cpp}
    paramHdlr.addAPIID(SomeAPI::getAPIIDString(),this) ;
    paramHdlr.addAPIID(ParamBEAPI::getAPIIDString(),&paramHdlr) ;
  \endcode

  Finally, load the parameters exported by SomeAPI into the handler.
  ParamBEAPI_Imp provides built-in support for parameters of simple type
  where the existing %API code provides dedicated get and set methods with
  signature
  \code{.cpp}
    ValType getParam () ;
    setParam(ValType val) ;
  \endcode
  Taking a parameter of type \p std::string as an example, the code to export
  the parameter is
  \code{.cpp}
    typedef ParamBEAPI_Imp<SomeAPI>::ParamEntry_Imp<std::string> SPE ;
    ParamBEAPI_Imp<SomeAPI>::ParamEntry *paramEntry ;
    paramEntry = new SPE("Param1",&SomeAPI::getParam1,&SomeAPI::setParam1) ;
    paramHdlr.addParam("Param1",paramEntry) ;
  \endcode

  \todo
  Provide built-in support for arbitrary parameter set/get methods by way of
  a signature
  \code{.cpp}
    void get(void *&blob)
    void set(void *&blob)
  \endcode
  The %API implementor would be responsible for writing a wrapper around %API
  methods to unpack / pack the blob.
*/
template<class Managed>
class ParamBEAPI_Imp : public ParamBEAPI {

public: 

  class ParamEntry ;

  /// ParamBEAPI_Imp implements ParamBEAPI.
  inline static const char *getAPIIDString ()
  { return (ParamBEAPI::getAPIIDString()) ; }

  /// \name Constructors and Destructor
  //@{
  /// Default constructor; \sa #create
  ParamBEAPI_Imp() ;
  /// Default copy constructor; \sa #clone
  ParamBEAPI_Imp(const ParamBEAPI_Imp &original) ;
  /// Virtual constructor
  ParamBEAPI *create() ;
  /// Virtual copy constructor
  ParamBEAPI *clone() ;
  /// Destructor
  ~ParamBEAPI_Imp() ;
  /// Assignment
  ParamBEAPI_Imp &operator=(const ParamBEAPI_Imp &rhs) ;
  //@}

/*! \name Back-end parameter reporting and manipulation methods.

  These methods are implemented by an API that supports parameter management.
  They provide the capability for the API to report the parameters that it
  exports and to get and set the values of those parameters.

  It's the responsibility of the API implementor to document the identifying
  strings for the exported parameters and the contents of blobs used to set
  and get values.
*/
//@{
  /*! \brief Return an array of identifiers for exported parameters.

    Each entry in the array should be a null-terminated string.
  */
  inline int reportParams (const char **&params)
  { params = paramIDs_ ;  return (paramIDCnt_) ; }

  /*! \brief Get the value of the requested parameter

    The value should be returned in the provided \p blob. A return value of
    false means failure to retrieve the value.
  */
  bool get(const char *ident, void *&blob) ;

  /*! \brief Set the value of the requested parameter

    Set the parameter to the value provided in \p blob. A return value of
    false means the parameter was not set.
  */
  bool set(const char *ident, const void *&blob) ;
//@}

/*! \name Utilities for use by the managed object

  Utility functions to allow the managed object to load information into this
  back-end parameter management object, and retrieve it.
*/
//@{
  /*! \brief Add a parameter to the list of exported parameters

    Note that ParamEntry is a base class for a set of templated subclasses. A
    pointer to \p entry is required here.
  */
  void addParam (const char *paramID, ParamEntry *entry) ;

  /*! \brief Add an API to the set of APIs implemented by the managed object

    Add the ident string for the %API and a pointer to the object that
    implements the %API. Where the managed object implements the %API
    by inheritance, \p obj will typically be \p this. Where the managed
    object implements the %API by composition, \p obj will typically be
    a member of the managed object.
  */
  void addAPIID (const char *apiID, void *obj) ;

  /*! \brief Retrieve the list of capability APIs stored by #addAPIID.

    Note that ParamBEAPI_Imp does not derive from API. This is not an
    override implementation for API::getAPIs (though it provides a trivial
    implementation for that method).
  */
  inline int getAPIs(const char **&idents)
  { idents = apiIDs_ ; return (apiIDCnt_) ; }

  /*! \brief Retrieve the object that implements a given %API.

    Note that ParamBEAPI_Imp does not derive from API. This is not an
    override implementation for API::getAPIPtr (though it provides a trivial
    implementation for that method).
  */
  inline void *getAPIPtr(const char *ident)
  { if (objForAPI_.find(ident) == objForAPI_.end())
    { return (nullptr) ; }
    else
    { return (objForAPI_.find(ident)->second) ; } }
//@}

/*! \name Exported parameter list entries

  The class ParamEntry and its templated derived classes ParamEntry_Imp
  record the information necessary to get and set the values of parameters in
  the managed object.

  Each instantiation of a templated class ParamEntry_Imp encodes the data
  type of its parameter; this makes it possible to write wrapper functions
  for get and set that properly handle the data blob. ParamEntry and its
  derived classes must be member classes of ParamBEAPI_Imp because the type
  of the managed object must be known in order to properly invoke the
  pointer-to-member-function pointers for the get and set methods.
*/
//@{
  /*! \brief Parameter list entry base class

    This base class provides a convenient common class for entries in the
    exported parameter list.
  */
  class ParamEntry
  { 
    public:
      virtual ParamEntry *clone() = 0 ;
      virtual ~ParamEntry() { }
    protected:
    /// Name of the parameter
    std::string paramID_ ;
    public:
    /// Interface to the get method
    virtual bool get (Managed *obj, void *&blob) = 0 ;
    /// Interface to the set method
    virtual bool set (Managed *obj, const void *&blob) = 0 ;
  } ;

  /*! \brief Parameter list entry

    An instance of this class is instantiated for a specific parameter data
    type (ValType) and a specific managed class (Managed). Together, this
    gives enough information to write get and set wrappers that will pack or
    unpack the data blob and invoke the proper method in the managed class
    by way of pointers-to-member-functions.
  */
  template<class ValType>
  class ParamEntry_Imp : public ParamEntry
  { 
    public:
    typedef ValType (Managed::*GetFunc)() const ;
    typedef void (Managed::*SetFunc)(ValType val) ;

    /// Initialising constructor
    ParamEntry_Imp (std::string id,
		    GetFunc getter, SetFunc setter)
    { ParamEntry::paramID_ = id ;
      getFunc_ = getter ;
      setFunc_ = setter ;
    } 
    /// Copy constructor
    ParamEntry_Imp (const ParamEntry_Imp &rhs)
    {
      ParamEntry::paramID_ = rhs.paramID_ ;
      getFunc_ = rhs.getFunc_ ;
      setFunc_ = rhs.setFunc_ ;
    }
    /// Virtual copy constructor
    ParamEntry *clone ()
    {
      return (new ParamEntry_Imp(*this)) ;
    }
    /// Destructor
    ~ParamEntry_Imp () { } ;

    /// Wrapper to invoke the get method from the managed class
    bool get (Managed *obj, void *&blob)
    {
      std::cout << "    invoking derived get." << std::endl ;
      ValType *val = static_cast<ValType *>(blob) ;
      *val = (obj->*getFunc_)() ;
      return (true) ;
    }

    /// Wrapper to invoke the set method from the managed class
    bool set (Managed *obj, const void *&blob)
    {
      std::cout << "    invoking derived set." << std::endl ;
      const ValType *val = static_cast<const ValType *>(blob) ;
      (obj->*setFunc_)(*val) ;
      return (true) ;
    }

    private:

    /// Pointer-to-member get function in the managed class
    GetFunc getFunc_ ;
    /// Pointer-to-member set function in the managed class
    SetFunc setFunc_ ;
  } ;
//@}


/* Begin private data for ParamBEAPI_Imp */

private:

  /// Actual number of capabilities
  int apiIDCnt_ ;
  /// Length of capability ID vector
  int apiIDLen_ ;
  /// Vector of capability API IDs
  const char **apiIDs_ ;

  /// API to implementing object map
  std::map<std::string,void*> objForAPI_ ;

  /// Actual number of parameters
  int paramIDCnt_ ;
  /// Length of parameter ID vector
  int paramIDLen_ ;
  /// Vector of parameter IDs
  const char **paramIDs_ ;

  typedef std::map<std::string,ParamEntry *> EntryForParamMap ;
  /// Param ID to entry map
  EntryForParamMap entryForParam_ ;

  /// Indicator; false if the message handler belongs to the client
  bool dfltHandler_ ;
  /// Message handler
  CoinMessageHandler *msgHandler_ ;
  /// Messages
  CoinMessages msgs_ ;
  /// Log (verbosity) level
  int logLvl_ ;

} ;

} // namespace Osi2

/*
  And then we need the code, so that we have the definitions handy when
  instantiating a particular instance of the templated classes.
*/
#include "Osi2ParamBEAPI_Imp.cpp"

#endif   // Osi2ParamBEAPI_Imp_HPP
