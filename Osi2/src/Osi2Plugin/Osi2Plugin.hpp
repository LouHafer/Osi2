/*
  Copyright 2011 Lou Hafer, Matt Saltzman
  This code is licensed under the terms of the Eclipse Public License (EPL)

  Based on original design and code by Gigi Sayfan published in five parts
  in Dr. Dobbs, starting November 2007.
*/
/*! \file Osi2Plugin.hpp
    \brief Type definitions for the Osi2 plugin manager.

  This file contains the type definitions (types, structures, enums) used by
  the plugin framework. Generally useful for the plugin manager and for
  plugins.
*/
#ifndef OSI2_PLUGIN_HPP
#define OSI2_PLUGIN_HPP

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

namespace Osi2 {

// Forward declarations
struct PlatformServices ;
struct RegisterParams ;
struct ObjectParams ;

/// A typedef for null-terminated strings, so they're easy to recognise
typedef uint8_t CharString ;

/// Languages (C, C++) supported by the plugin framework.
enum PluginLang
{ Plugin_C, Plugin_CPP } ;

/// Plugin manager API version number
struct PluginAPIVersion
{
  int32_t major_ ;
  int32_t minor_ ;
} ;

/*! \brief A plugin state object

  This object is responsible for keeping track of whatever information is
  necessary to manage the plugin's state, if any. The plugin manager allows
  (but does not require) a single state object for the plugin library and
  individual state objects for each API. These can all be different objects
  or one and the same object; whatever is appropriate for the plugin. As
  far as the manager is concerned, they are opaque pointers.

  If the plugin library requires a library state object, a pointer to
  this object should be stored in the PlatformServices object passed as a
  parameter to the \link Osi2::InitFunc initialisation function \endlink.
  The pointer will be handed back to the plugin in the PlatformServices
  object provided as part of the parameters of any call made by the
  PluginManager to the plugin (the functions \link Osi2::CreateFunc \endlink,
  \link Osi2::DestroyFunc \endlink, and \link Osi2::ExitFunc \endlink.

  If an API requires an API state object, a pointer to this object should
  be stored in the RegisterParams object passed as a parameter to the
  \link Osi2::RegisterFunc registration function \endlink. The pointer
  will be handed back to the plugin in the ObjectParams object provided
  as part of the parameters to the object \link Osi2::CreateFunc creation
  \endlink and \link Osi2::DestroyFunc deletion \endlink functions.

  In a typical use case, the plugin state object will be an object of the
  class that defines the (static) functions used for communication between
  the plugin manager and the plugin. In essence, we're handing these static
  functions an explicit pointer to an object of their defining class which
  can be used to maintain state information.
*/
typedef void PluginState ;

/*! \brief Unique ID for a plugin library

  In order to match API registrations to plugin libraries, a unique ID is
  needed. This ID is generated by the plugin manager when the plugin library
  is loaded. It is passed to the plugin library in the PlatformServices
  parameter of the \link Osi2::InitFunc initialisation function \endlink.
  It must be returned by the plugin library in the RegisterParam block
  when the plugin \link Osi2::RegisterFunc registers an API \endlink. This
  allows the manager to unambiguously match an API registration to the plugin
  library that provides the API.

  In the current implementation, this is just a reference to the
  DynamicLibrary object for the library.
*/
typedef void *PluginUniqueID ;

/*! \defgroup PluginFuncs Plugin Management Functions

  These functions are used for communication between the plugin and the
  plugin manager.
  
  The functions \link Osi2::RegisterFunc \endlink and \link
  Osi2::InvokeServiceFunc \endlink are defined by
  the PluginManager and passed to the plugin in the PlatformServices
  parameter of the \link Osi2::InitFunc initialisation function \endlink.
  
  All other functions are defined by the plugin library. They must be defined
  by any plugin library, even if they do nothing.
*/
//@{
/*! \brief Constructor for a plugin object

  This function is implemented by the plugin and invoked by the PluginManager
  to create instances of objects (APIs) supported by the plugin. A
  \c CreateFunc must be specified for each registered API, but exactly
  how the function is implemented is a decision for the plugin author.

  \param parms \link ObjectParams object parameters \endlink for use
  	   by the plugin.
  \returns Pointer to an object supporting the requested API, or null for
	   failure.
*/
typedef void *(*CreateFunc)(const ObjectParams *parms) ;

/*! \brief Destructor for a plugin object

  This function is implemented by the plugin and invoked by the PluginManager
  to destroy instances of objects (APIs) supported by the plugin. A
  \c DestroyFunc must be specified for each registered API, but exactly how
  the function is implemented is a decision for the plugin author.

  \param victim the object to be destroyed.
  \param parms \link Osi2::ObjectParams object parameters \endlink for use
  	   by the plugin.
  \returns 0 for success, nonzero for failure.
*/
typedef int32_t (*DestroyFunc)(void *victim, const ObjectParams *parms) ;

/*! \brief API registration function

  This function is implemented by the PluginManager and passed to the plugin
  in the PlatformServices parameter of the \link Osi2::InitFunc
  initialisation function \endlink. The plugin should invoke this function to
  register each API that it can provide. Typically this method is called as
  part of the execution of the \link Osi2::InitFunc initialisation function
  \endlink. It can also be invoked at any time while the plugin is active to
  register new APIs.

  \param apiStr A null-terminated character string identifying the API.
  \param parms  \link Osi2::RegisterParams registration information \endlink
		 for consumption by the plugin manager.
  \returns 0 for a successful registration, nonzero otherwise.
*/
typedef int32_t (*RegisterFunc)(const CharString *apiStr,
				const RegisterParams *parms) ;

/*! \brief Function to allow the plugin to invoke services provided by the
	   plugin manager.

  This method is implemented by the PluginManager and passed to the plugin in
  a PlatformServices parameter object. The plugin can invoke this function to
  access services provided by the PluginManager.

  \param serviceName A null-terminated character string identifying the
		     service.
  \param serviceParams arbitrary parameter block appropriate to service

  \returns 0 if the service is successfully invoked, nonzero otherwise

  \todo Implement this! Currently unimplemented.
*/
typedef int32_t (*InvokeServiceFunc)(const CharString *serviceName,
					  void *serviceParams) ;

/*! \brief Type definition of the \c exitPlugin function

  This function is called by the PluginManager to tell the plugin to clean
  up in preparation for unloading the plugin library. The plugin library's
  control object (as returned by the \link Osi2::InitFunc initialisation
  function \endlink will be included, as will the library's \link
  PluginUniqueID unique id \endlink.

  The PluginManager will locate the \c exitPlugin function with
  DynamicLibrary::getSymbol using the string \c "exitPlugin" during plugin
  initialisation.

  \param services Services provided by the plugin manager.
  \returns 0 for successful execution, nonzero otherwise.
*/
typedef int32_t (*ExitFunc)(const PlatformServices *services) ;

/*! \brief Type definition of the \c initPlugin function

  The plugin manager will locate the \c initPlugin function with
  DynamicLibrary::getSymbol using the string \c "initPlugin" and call it
  to initialise the plugin library.

  If the plugin library requires a state object, the initialisation function
  should fill in the \link Osi2::PlatformServices#ctrlObj_ ctrlObj_ \endlink
  field in the \p services parameter.

  If initialisation fails for any reason, the plugin is expected to take any
  actions required for clean up (in essence, ready itself for library unload)
  and return null. It may also (optionally) report the error via the error
  reporting function provided in the PlatformServices parameter.

  \param services Services provided by the plugin manager.
  \returns A pointer to the library's \link Osi2::ExitFunc exit function
	   \endlink, or null if initialisation fails.
 */
typedef ExitFunc (*InitFunc)(const PlatformServices *services) ;

//@}


/*! \brief Parameters required to register an API with the plugin manager

  The plugin supplies this structure when it registers with the plugin
  manager.
*/
struct RegisterParams
{
  /*! \name Management information

    This information is used by the plugin manager for validation and
    management.
  */
  /// Plugin manager version expected by the plugin
  PluginAPIVersion version_ ;

  /*! \brief Unique plugin ID

    Supplied by the plugin manager in the \link Osi2::PlatformServices
    parameters \endlink for the \link Osi2::InitFunc initialisation function
    \endlink when the plugin library is initialised.  Must be specified when
    registering an API in order to associate the API to the plugin library.
  */
  PluginUniqueID pluginID_ ;

  /// Implementation language of the API
  PluginLang lang_ ;

  /*! \brief Plugin control object for this API

    Optional pointer to an object holding state information used to manage
    objects supporting this API. This reference will be returned to the plugin
    as part of the \link Osi2::ObjectParams parameters \endlink to the
    \link Osi2::CreateFunc create \endlink and \link Osi2::DestroyFunc destroy
    \endlink functions for objects of this API.
  */
  PluginState *ctrlObj_ ;

  /// Constructor for API being registered
  CreateFunc createFunc_ ;
  /// Destructor for API being registered
  DestroyFunc destroyFunc_ ;
} ;

/*! \brief Parameters passed to a plugin for object creation and destruction.

  The functions to create and destroy an object are static methods, but for
  plugins capable of providing multiple APIs it may be useful to pass a
  pointer to the plugin object.
*/
struct ObjectParams
{
  /// Character string specifying the API(s) to be supplied
  const CharString *apiStr_ ;
  /*! \brief Services provided by the plugin manager

    See Osi2::PlatformServices.
  */
  const struct PlatformServices *platformServices_ ;
  /*! \brief Plugin control object for this API

    Optional pointer to an object holding state information used to manage
    objects supporting this API. This is the control object supplied with the
    \link Osi2::RegisterParams registration parameters \endlink for this API.
  */
  PluginState *ctrlObj_ ;
} ;

/*! \brief Attributes and methods provided by the plugin framework for use
	   by plugins.

  A reference to this structure is provided, directly or indirectly, in the
  parameters of all functions used for communication between the plugin
  manager and the plugin library. In general, there are two classes of
  information: information about the plugin manager, useable by all plugin
  libraries; and information specific to the plugin library (the library's
  unique ID and control state object).

  \todo Provide some services!
*/
struct PlatformServices
{
  /// Plugin manager version
  PluginAPIVersion version_ ;
  /// Default directory to search for plugins
  const CharString *dfltPluginDir_ ;
  /// Unique ID assigned to this plugin library
  PluginUniqueID pluginID_ ;
  /// Control object for this plugin library
  PluginState *ctrlObj_ ;
  /// Method to register an object supported by the plugin
  RegisterFunc registerObject_ ; 
  /*! \brief Method to allow the plugin to invoke a service provided by the
  	     plugin manager.
  */
  InvokeServiceFunc invokeService_; 
} ;


}  // end namespace Osi2

#ifdef  __cplusplus
}
#endif

#endif

