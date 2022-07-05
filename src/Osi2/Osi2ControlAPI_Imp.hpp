/*
  Copyright 2011 Lou Hafer, Matt Saltzman
  This code is licensed under the terms of the Eclipse Public License (EPL)

  $Id$
*/
/*! \file Osi2ControlAPI_Imp.hpp

  Defines an implementation of Osi2::ControlAPI, an API to manage
  plugin libraries and objects.
*/


#ifndef Osi2ControlAPI_Imp_HPP
#define Osi2ControlAPI_Imp_HPP

#include "Osi2PluginManager.hpp"

#include "Osi2ParamBEAPI_Imp.hpp"
#include "Osi2ControlAPI.hpp"
#include "Osi2CtrlAPIMessages.hpp"

namespace Osi2 {

/*! \brief Osi2 ControlAPI implementation

  This class implements the Osi2 Control API (Osi2::ControlAPI).
*/
class OSI2LIB_EXPORT ControlAPI_Imp : public ControlAPI {

public:

  /// ControlAPI_Imp implements ControlAPI.
  inline static const char *getAPIIDString ()
  { return (ControlAPI::getAPIIDString()) ; }

    /// \name Constructors and Destructor
    //@{
    /// Virtual constructor
    ControlAPI *create() ;
    /// Virtual copy constructor
    ControlAPI *clone() ;
    /// Destructor
    ~ControlAPI_Imp() ;
    /// Default constructor; \sa #create
    ControlAPI_Imp() ;
    /// Default copy constructor; \sa #clone
    ControlAPI_Imp(const ControlAPI_Imp &original) ;
    /// Assignment
    ControlAPI_Imp &operator=(const ControlAPI_Imp &rhs) ;
    //@}

    /*! \name Library load and unload
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
    virtual int load(const std::string &shortName) ;

    /*! \brief Load the specified plugin library

      Attempt to load library \p libName from the default plugin library
      directory. Associate the plugin library with \p shortName.
    */
    virtual int load(const std::string &shortName, const std::string &libName) ;

    /*! \brief Load the specified plugin library

      Attempt to load library \p libName from directory \p dirName. Associate
      the plugin library with \p shortName.
    */
    virtual int load(const std::string &shortName,
                     const std::string &libName, const std::string *dirName) ;

    /*! \brief Unload the specified library.

      Unloads the specified library. The return value will be 0 if all went
      smoothly, non-zero on error.
    */
    virtual int unload(const std::string &shortName) ;

    /// Find the short name for the specified library
    virtual std::string getShortName(PluginUniqueID libID) ;

    /// Find the full path for the specified library
    virtual std::string getFullPath(PluginUniqueID libID) ;
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
      API. If \p shortName is specified but not registered, the method will
      issue a warning and proceed as if no restriction was given.

      If the call executes without error, \p obj will contain a reference to an
      object supporting the requested API, otherwise it will be set to null.

      \returns:
        -1: error
         0: the object was successfully created.
         1: the object was successfully created but the plugin restriction was
      ignored because it was not recognised

    */
    virtual int createObject(API *&obj, const std::string &apiName,
                             const std::string *shortName = 0) ;

    /*! \brief Destroy the specified object

      In general, invoking delete on the object will work just fine. This method
      is necessary only if the plugin providing the object needs to know about
      its demise.

      \p obj will be set to null on return.

      \returns:
        -1: an error occurred during destruction
         0: destruction completed without error
    */
    virtual int destroyObject(API *&obj) ;

    //@}

    /*! \name Control API control methods

      Miscellaneous methods that control the behaviour of a ControlAPI object.
    */
    //@{

    /// Get the default plugin directory
    inline std::string getDfltPluginDir() const {
        return (dfltPluginDir_) ;
    }

    /// Set the default plugin directory
    inline void setDfltPluginDir(std::string dfltDir) {
        dfltPluginDir_ = dfltDir ;
    }


    /// Set the log (verbosity) level
    inline void setLogLvl(int logLvl) {
        logLvl_ = logLvl ;
        msgHandler_->setLogLevel(logLvl_) ;
    }

    /// Get the log (verbosity) level
    inline int getLogLvl() const {
        return (logLvl_) ;
    }

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

    /// Parameter management forwarding functions
    inline int getAPIs(const char **&idents)
    { return (paramHandler_.getAPIs(idents)) ; }
    inline void *getAPIPtr (const char *ident)
    { return (paramHandler_.getAPIPtr(ident)) ; }

    //@}

private:

    /// Cached reference to plugin manager.
    PluginManager *pluginMgr_ ;

    /*! \brief Plugin library management information

      This data structure holds the information that the ControlAPI uses to
      manage plugin libraries.
    */
    struct DynLibInfo {
        /// Full path for the library
        std::string fullPath_ ;
        /// Unique ID for the library
        PluginUniqueID uniqueID_ ;
    } ;
    /// Map type for knownLibMap_
    typedef std::map<std::string, DynLibInfo> LibMapType ;
    /// Map to associate short names with full paths for plugin libraries
    LibMapType knownLibMap_ ;

    /// Default plugin library directory
    std::string dfltPluginDir_ ;

    /// Indicator; false if the message handler belongs to the client
    bool dfltHandler_ ;
    /// Message handler
    CoinMessageHandler *msgHandler_ ;
    /// Messages
    CoinMessages msgs_ ;
    /// Log (verbosity) level
    int logLvl_ ;

    /// Parameter management object
    ParamBEAPI_Imp<ControlAPI> paramHandler_ ;

    /*! \brief Control information

      This class defines how API object control information is structured for
      this implementation of the control API.
    */
    struct APIObjCtrlInfo {
      /// Initialising constructor
      APIObjCtrlInfo(std::string apiName, PluginUniqueID libID)
        : apiName_(apiName),
	  libID_(libID)
      {}

      /// API name
      const std::string apiName_ ;
      /// Library unique ID
      const PluginUniqueID libID_ ;
    } ;

} ;

} // namespace Osi2 ;

#endif
