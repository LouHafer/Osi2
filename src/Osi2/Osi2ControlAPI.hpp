
/*
  Copyright 2011 Lou Hafer, Matt Saltzman
  This code is licensed under the terms of the Eclipse Public License (EPL)
*/

/*! \file Osi2ControlAPI.hpp

  Defines the virtual class Osi2::ControlAPI, an API to manage plugin
  libraries and objects.
*/


#ifndef Osi2ControlAPI_HPP
# define Osi2ControlAPI_HPP

#include <string>

#include "Osi2API.hpp"

namespace Osi2 {

/*! \brief Osi2 ControlAPI virtual base class

  This abstract class defines the Osi2 control interface, including loading and
  unloading of objects.
*/
class ControlAPI : public API {

public:

    /// Return the "ident" string for the control %API.
    inline static const char *getAPIIDString () { return ("ControlAPI") ; }

    /// \name Constructors and Destructors
    //@{
    /// Virtual constructor
    virtual ControlAPI *create() = 0 ;
    /// Virtual copy constructor
    virtual ControlAPI *clone() = 0 ;
    /// Virtual destructor
    virtual ~ControlAPI() {}
    //@}

    /*! \name Library Management
        \brief Methods to load and unload plugin libraries

      Given only a short name xxx the method will look for libOsi2XxxShim.so in
      the default plugin library directory. The other load methods allow for an
      arbitrary mapping between the short name and the plugin library.
    */
    //@{

    /*! \brief Load the specified plugin library

      Given \p shortName 'xxx', attempt to load libOsi2XxxShim.so from the
      default plugin library directory. Associate the plugin library with
      \p shortName.
    */
    virtual int load(const std::string &shortName) = 0 ;

    /*! \brief Load the specified plugin library

      Attempt to load library \p libName from the default plugin library
      directory. Associate the plugin library with \p shortName.
    */
    virtual int load(const std::string &shortName,
                     const std::string &libName) = 0 ;

    /*! \brief Load the specified plugin library

      Attempt to load library \p libName from directory \p dirName. Associate
      the plugin library with \p shortName.
    */
    virtual int load(const std::string &shortName,
                     const std::string &libName, const std::string *dirName) = 0 ;

    /*! \brief Unload the specified library.

      Unloads the specified library. The return value will be 0 if all went
      smoothly, non-zero on error.
    */
    virtual int unload(const std::string &shortName) = 0 ;

    //@}

    /*! \name API (Object) Management
        \brief Methods to create and destroy API objects.
    */
    //@{

    /*! \brief Create an object of the specified API

      The parameter \p apiName should contain a string specifying the API. This
      string must be known to some plugin library. (In general, the APIs
      supported by a plugin will be documented with the plugin.)

      If \p shortName is specified, only the specified plugin library will be
      considered. If \p shortName is not specified, an arbitrary choice will be
      made from plugin libraries capable of supplying an object supporting the
      API.

      If the call executes without error, \p obj will contain a reference to an
      object supporting the requested API, otherwise it will be set to null.

      \returns:
        -1: error
         0: the object was successfully created.

    */
    virtual int createObject(API *&obj, const std::string &apiName,
                             const std::string *shortName = 0) = 0 ;

    /*! \brief Destroy the specified object

      In general, invoking delete on the object will work just fine. This method
      is necessary only if the plugin providing the object needs to know about
      its demise.
    */
    virtual int destroyObject(API *&obj) = 0 ;

    //@}

    /*! \name Control API control methods

      Miscellaneous methods that control the behaviour of a ControlAPI object.
    */
    //@{

    /// Get the default plugin directory
    virtual std::string getDfltPluginDir() const = 0 ;

    /// Set the default plugin directory
    virtual void setDfltPluginDir(std::string dfltDir) = 0 ;

    /// Set the log (verbosity) level
    virtual void setLogLvl(int logLvl) = 0 ;

    /// Get the log (verbosity) level
    virtual int getLogLvl() const = 0 ;

    //@}

} ;

} // namespace Osi2 ;

#endif
