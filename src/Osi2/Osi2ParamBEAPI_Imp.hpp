/*
  Copyright 2019 Lou Hafer, Matt Saltzman, Kevin Poton
  This code is licensed under the terms of the Eclipse Public License (EPL)

  $Id$
*/
/*! \file Osi2ParamBEAPI_Imp.hpp

  Defines Osi2::ParamBEAPI_Imp, an implementation of the Osi2::ParamBEAPI
  interface.
*/


#ifndef Osi2ParamBEAPI_Imp_HPP
#define Osi2ParamBEAPI_Imp_HPP

#include <map>

#include "Osi2API.hpp"
#include "Osi2APIMgmt_Imp.hpp"

#include "Osi2ParamBEAPI.hpp"
#include "Osi2ParamBEAPIMessages.hpp"

namespace Osi2 {

/*! \brief An implementation of the parameter management back end.

  This class implements the ParamBEAPI methods that must be supplied by
  an %API implementation to support parameter management. It is intended
  to be used by including a ParamBEAPI_Imp object as a member in the client
  %API object.

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
  This isn't integral to the handler's operation; it's a convenience service
  provided to the client. (The client must support %API inquiry via
  derivation from Osi2::API; this service makes that trivial. \sa
  Osi2::APIMgmt.)

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
template<class Client>
class ParamBEAPI_Imp : public ParamBEAPI {

public: 

  class ParamEntry ;

  /// ParamBEAPI_Imp implements ParamBEAPI.
  inline static const char *getAPIIDString ()
  { return (ParamBEAPI::getAPIIDString()) ; }

  /// \name Constructors and Destructor
  //@{
  /// constructor; \sa #create
  ParamBEAPI_Imp(Client *client) ;
  /// Default copy constructor; \sa #clone
  ParamBEAPI_Imp(const ParamBEAPI_Imp &original) ;
  /// Virtual constructor
  ParamBEAPI *create(Client *client) ;
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

/*! \name Utilities for use by the client object

  Utility functions to allow the client object to load information into this
  back-end parameter management object, and retrieve it.
*/
//@{
  /*! \brief Add a parameter to the list of exported parameters

    Note that ParamEntry is a base class for a set of templated subclasses. A
    pointer to \p entry is required here.
  */
  void addParam (const char *paramID, ParamEntry *entry) ;

  /*! \brief Add an API to the set of APIs implemented by the client object

    Add the ident string for the %API and a pointer to the object that
    implements the %API. Where the client object implements the %API
    by inheritance, \p obj will typically be \p this. Where the client
    object implements the %API by composition, \p obj will typically be
    a member of the client object.

    Provided as a convenience; \sa Osi2::APIMgmt_Imp.
  */
  inline void addAPIID (const char *apiID, void *obj)
  { apiMgr_.addAPIID(apiID,obj) ; }

  /*! \brief Retrieve the list of capability APIs stored by #addAPIID.

    Note that ParamBEAPI_Imp does not derive from API. This is not an
    override implementation for API::getAPIs (though it provides a trivial
    implementation for that method).

    Provided as a convenience; \sa Osi2::APIMgmt_Imp.
  */
  inline int getAPIs(const char **&idents)
  { return (apiMgr_.getAPIs(idents)) ; }

  /*! \brief Retrieve the object that implements a given %API.

    Note that ParamBEAPI_Imp does not derive from API. This is not an
    override implementation for API::getAPIPtr (though it provides a trivial
    implementation for that method).

    Provided as a convenience; \sa Osi2::APIMgmt_Imp.
  */
  inline void *getAPIPtr(const char *ident)
  { return (apiMgr_.getAPIPtr(ident)) ; }
//@}

/*! \name Exported parameter list entries

  The class ParamEntry and its templated derived classes ParamEntry_Imp
  and ParamEntry_Gen record the information necessary to get and set the
  values of parameters in the client object.

  Each instantiation of one of the templated classes encodes the data
  type of its parameter; this makes it possible to write wrapper functions
  for get and set that properly handle the data blob. ParamEntry and its
  derived classes must be member classes of ParamBEAPI_Imp because the type
  of the client object must be known in order to properly invoke the
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
    virtual bool get (Client *obj, void *&blob) = 0 ;
    /// Interface to the set method
    virtual bool set (Client *obj, const void *&blob) = 0 ;
  } ;

  /*! \brief Parameter list entry (specific get/set)

    This class handles parameters with specific get & set methods (i.e., a
    dedicated pair of methods for each parameter).

    An instance of this class is instantiated for a specific parameter data
    type (ValType) and a specific client class (Client). Together, this
    gives enough information to write get and set wrappers that will pack or
    unpack the data blob and invoke the proper method in the client class
    by way of pointers-to-member-functions.
  */
  template<class ValType>
  class ParamEntry_Imp : public ParamEntry
  { 
    public:
    typedef ValType (Client::*GetFunc)() const ;
    typedef void (Client::*SetFunc)(ValType val) ;

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

    /// Wrapper to invoke the get method from the client class
    bool get (Client *obj, void *&blob)
    {
      ValType *val = static_cast<ValType *>(blob) ;
      *val = (obj->*getFunc_)() ;
      return (true) ;
    }

    /// Wrapper to invoke the set method from the client class
    bool set (Client *obj, const void *&blob)
    {
      const ValType *val = static_cast<const ValType *>(blob) ;
      (obj->*setFunc_)(*val) ;
      return (true) ;
    }

    private:

    /// Pointer-to-member get function in the client class
    GetFunc getFunc_ ;
    /// Pointer-to-member set function in the client class
    SetFunc setFunc_ ;
  } ;

  /*! \brief Parameter list entry (generic get/set)

    This class handles parameters that have generic get and set methods
    (i.e., a pair of methods that take an identifier of some sort to
    specify the parameter).

    An instance of this class is instantiated for a specific parameter data
    type (ValType) and a specific client class (Client). Together, this
    gives enough information to write get and set wrappers that will pack or
    unpack the data blob and invoke the proper method in the client class
    by way of pointers-to-member-functions.
  */
  template<class ValType>
  class ParamEntry_Gen : public ParamEntry
  { 
    public:
    typedef ValType (Client::*GetFunc)(std::string) const ;
    typedef void (Client::*SetFunc)(std::string, ValType val) ;

    /// Initialising constructor
    ParamEntry_Gen (std::string id,std::string what,
		    GetFunc getter, SetFunc setter)
    { ParamEntry::paramID_ = id ;
      what_ = what ;
      getFunc_ = getter ;
      setFunc_ = setter ;
    } 
    /// Copy constructor
    ParamEntry_Gen (const ParamEntry_Gen &rhs)
    {
      ParamEntry::paramID_ = rhs.paramID_ ;
      what_ = rhs.what_ ;
      getFunc_ = rhs.getFunc_ ;
      setFunc_ = rhs.setFunc_ ;
    }
    /// Virtual copy constructor
    ParamEntry *clone ()
    {
      return (new ParamEntry_Gen(*this)) ;
    }
    /// Destructor
    ~ParamEntry_Gen () { } ;

    /// Wrapper to invoke the get method from the client class
    bool get (Client *obj, void *&blob)
    {
      ValType *val = static_cast<ValType *>(blob) ;
      *val = (obj->*getFunc_)(what_) ;
      return (true) ;
    }

    /// Wrapper to invoke the set method from the client class
    bool set (Client *obj, const void *&blob)
    {
      const ValType *val = static_cast<const ValType *>(blob) ;
      (obj->*setFunc_)(what_,*val) ;
      return (true) ;
    }

    private:

    /// Parameter specifier
    std::string what_ ;
    /// Pointer-to-member get function in the client class
    GetFunc getFunc_ ;
    /// Pointer-to-member set function in the client class
    SetFunc setFunc_ ;
  } ;

  /*! \brief Parameter list entry (pass through)

    This class handles parameters that use complex parameter structures. It
    assumes that the user will supply a pointer to a blob of space which
    is passed through to the set / get methods of the client class.

    An instance of this class is instantiated for a specific client class
    (Client); the data type is handled as void*.  This provides enough
    information to invoke client methods. It's up to the client's methods
    to know how to unpack or pack the blob.
  */
  class ParamEntry_Void : public ParamEntry
  { 
    public:
    typedef void (Client::*GetFunc)(void *) const ;
    typedef void (Client::*SetFunc)(const void *) ;

    /// Initialising constructor
    ParamEntry_Void (std::string id,
		    GetFunc getter, SetFunc setter)
    { ParamEntry::paramID_ = id ;
      getFunc_ = getter ;
      setFunc_ = setter ;
    } 
    /// Copy constructor
    ParamEntry_Void (const ParamEntry_Void &rhs)
    {
      ParamEntry::paramID_ = rhs.paramID_ ;
      getFunc_ = rhs.getFunc_ ;
      setFunc_ = rhs.setFunc_ ;
    }
    /// Virtual copy constructor
    ParamEntry *clone ()
    {
      return (new ParamEntry_Void(*this)) ;
    }
    /// Destructor
    ~ParamEntry_Void () { } ;

    /// Wrapper to invoke the get method from the client class
    bool get (Client *obj, void *&blob)
    {
      (obj->*getFunc_)(blob) ;
      return (true) ;
    }

    /// Wrapper to invoke the set method from the client class
    bool set (Client *obj, const void *&blob)
    {
      (obj->*setFunc_)(blob) ;
      return (true) ;
    }

    private:

    /// Pointer-to-member get function in the client class
    GetFunc getFunc_ ;
    /// Pointer-to-member set function in the client class
    SetFunc setFunc_ ;
  } ;

//@}


/* Begin private data for ParamBEAPI_Imp */

private:

  /// The object we're serving
  Client *client_ ;
  /// API management object
  APIMgmt_Imp apiMgr_ ;

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
