#ifndef OSI2_DYNAMIC_LIBRARY_HPP
#define OSI2_DYNAMIC_LIBRARY_HPP

/*
  Copyright 2011 Lou Hafer, Matt Saltzman
  This code is licensed under the terms of the Eclipse Public License (EPL)

  Based on original design and code by Gigi Sayfan published in five parts
  in Dr. Dobbs, starting November 2007.
*/
/*! \file Osi2DynamicLibrary.hpp

  Declarations for the underlying dynamic library and symbol management
  methods.
*/

#include <string>
#include "Osi2Config.h"

namespace Osi2 {

/*! \brief A class to manage dynamic library load and unload.

  A call to the static method #load will load a library and return a
  DynamicLibrary object. The #getSymbol method can then be used to load
  individual methods from the library. Unloading the library is handled by
  the destructor for the DynamicLibrary object.

  At present, this library supports only dlopen / dlsym / dlclose (dlfcn.h).

  \todo Restore/check Windows functionality.

  \todo Convert/augment with libltdl.

  \todo Should there be a distinct unload method to invoke dlclose? Currently
        the destructor is the only way to close.

  \todo Consider making this a virtual class, with implementation classes for
        vanilla dlopen/dlclose, libtool ltdl, Windows, etc.
*/
/*
  Doxygen's builtin notion of grouping gets in the way here. What you see is a
  compromise.  -- lh, 110912 --
*/
class OSI2LIB_EXPORT DynamicLibrary {

public:

/*! \name Library Management
    \brief Library management methods

    Utilities to load a library. Unloading the library is handled via the
    destructor for the DynamicLibrary object.
*/
//@{

  /*! \brief Load a plugin

    This method will load the plugin specified in \p path. If successful,
    it will return a pointer to a DynamicLibrary object for the plugin. If
    unsuccessful, it will return nullptr and \p errStr will be loaded
    with an error message.
  */
  static DynamicLibrary *load(const std::string &path, std::string &errStr) ;
//@}

  /*! \name Symbol Management
      \brief Symbol management methods
  */
//@{
  /*! \brief Load a symbol with type conversion

    Load the specified symbol with appropriate type, using templating to
    avoid an explicit cast.  If unsuccessful, it will return nullptr and
    \p errStr will be loaded with an error message.
  */
  template<class T>
  T getSymbol(const std::string &name, std::string &errStr) ;

  /* \brief Load a symbol

    Load the specified symbol. The return value must be cast to the
    appropriate type for use. If unsuccessful, it will return nullptr and
    \p errStr will be loaded with an error message.

    Mainly here to avoid breaking existing code. May go away at some point.
    --lh, 190804 --
  */
  inline void *getSymbol(const std::string &name, std::string &errStr)
  { return (getSymbol<void *>(name,errStr)) ; }

  /*! \brief Load a function symbol

    Converting the pointer-to-object (void *) returned by dlsym to a
    pointer-to-function with arbitrary function signature requires beating
    C++ about the head with a 2x4. This method hides that ugliness. If
    unsuccessful, it will return nullptr and \p errStr will be loaded with
    an error message.
   */
  template<class FuncSig>
  FuncSig getFunc(const std::string &name, std::string &errStr) ;

//@}

    /*! \name Miscellaneous
    */
//@{
    /// Get the full path for the library
    inline std::string getLibPath () const {
        return (fullPath_) ;
    }
    /// Convenient stream i/o operator
    friend
    std::ostream &operator<<(std::ostream &ostr, const DynamicLibrary *lib) ;
//@}

    /*! \name Destructor */
//@{
    /*! \brief Destructor

      Unloads the library as part of execution.
    */
    ~DynamicLibrary() ;
//@}

private:

    /*! \name Constructors
        \brief Constructors

      Constructors are private because there should be exactly one instance of a
      DynamicLibrary object associated with a successfully loaded library. This
      object is returned by #load.
    */
//@{
    /// Default constructor
    DynamicLibrary() ;

    /*! \brief Constructor from dynamic library handle

      The parameter is the platform-specific dynamic library handle.
    */
    DynamicLibrary(void *handle) ;

    /// Copy constructor
    DynamicLibrary(const DynamicLibrary &dynlib) ;
//@}

    /// Platform-specific dynamic library handle.
    void *handle_;

    /// Full path of associated library
    std::string fullPath_ ;

    /// Default plugin directory
    std::string dfltPluginDir_ ;
} ;

}  // end namespace Osi2

/*
  The body for the templated methods has to be around at the point of
  instantiation.
*/

#include <iostream>
#include <dlfcn.h>

namespace Osi2 {

template<class T>
T DynamicLibrary::getSymbol (const std::string &name,
			     std::string &errorString)
{
  if (handle_ == nullptr) return (nullptr) ;

  void *sym = ::dlsym(handle_,name.c_str()) ;
  if (sym == nullptr) {
    errorString += "Failed to load symbol \"" + name + '"' ;
    const char *zErrorString = ::dlerror() ;
    if (zErrorString) errorString = errorString + ": " + zErrorString ;
  }
  return (static_cast<T>(sym)) ;
}

/*
  FuncType func =
    reinterpret_cast<FuncType>(lib->getSymbol(funcName,errStr)) ;

  The code above is what should really happen here. But up through C++11,
  ISO C++ required a warning for conversion of pointer-to-object to
  pointer-to-function. There are systems out there that may not support
  this conversion.  POSIX compatible systems must support conversion of
  pointer-to-object to pointer-to-function, otherwise dlsym (which returns
  a void*) simply won't work.

  C++ guarantees we can cast pointer-to-anything to size_t and get it back
  (with the possible exception of systems with segmented address space).
  So we slip through a blind spot in the compiler's algorithm for generating
  warnings.
*/
template <class FuncSig>
FuncSig DynamicLibrary::getFunc (const std::string &name,
				 std::string &errStr)
{
  size_t grossHack = 
      reinterpret_cast<size_t>(getSymbol(name,errStr)) ;
  FuncSig func = reinterpret_cast<FuncSig>(grossHack) ;
  return (func) ;
}

}  // end namespace Osi2

#endif
