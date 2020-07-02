/*
  Copyright 2011 Lou Hafer, Matt Saltzman
  This code is licensed under the terms of the Eclipse Public License (EPL)

  $Id$
*/
/*! \file Osi2CFuncTemplates.cpp
    \brief Templated functions for invoking dynamically loaded C functions.

  In the templates,
    - ObjType is the type of the object being acted on.
    - RetType is the return type of the C function
    - ValType is the type of the parameter to the C function

  To invoke the templated methods, the parameters are
    - lib: a DynamicLibrary object
    - obj: the object being acted on
    - funcName: the name of the function to be invoked
*/

/*
  Templates for simple get and set methods, templated for the object to be
  acted on and the type of the parameter. The heart of the template is the
  code to dynamically load the function, if required, and invoke it on the
  object. All the various overloads cope with passing parameters to the
  function to be invoked.

  The function signatures FuncSig are as follows:

  Simple getter, no parameters other than the object to be acted on:
    RetType (*FuncSig)(ObjType*)
  Simple getter, up to two additional parameters
    RetType (*FuncSig)(ObjType*,ValType1)
    RetType (*FuncSig)(ObjType*,ValType1,Valtype2)

  Simple setter, up to three parameters to pass in a value:
    void (*FuncSig)(ObjType *)
    void (*FuncSig)(ObjType *,ValType1)
    void (*FuncSig)(ObjType *,ValType1,ValType2)
    void (*FuncSig)(ObjType *,ValType1,ValType2,ValType3)
  The signature with no ValType is a simple imperative (e.g., doAction). The
  others allow for parameters (e.g., setValue(parm1))

  Don't get hung up on the simpleGetter / simpleSetter naming convention.
  Really, the choice comes down to whether or not the underlying C method
  returns a value.
  

  \todo Think of a way to implement a cache for function addresses so we're
  not invoking dlsym every time we use a function. Really, this cache ought to
  be a service of DynamicLibrary. Keep a map in the DynamicLibrary object?
*/

using Osi2::DynamicLibrary ;

template <class ObjType, class RetType>
RetType simpleGetter (DynamicLibrary *lib, ObjType *obj,
	      std::string funcName)
{ std::string errStr ;
  typedef RetType (*FuncSig)(ObjType *) ;
  FuncSig getterFunc =
    lib->getFunc<FuncSig>(funcName,errStr) ;
  if (getterFunc != nullptr) {
    return (getterFunc(obj)) ;
  } else {
    std::cout
      << "  simpleGetter(0): lookup failed for \"" << funcName << "\": "
      << errStr << std::endl ;
  }
  return (RetType()) ;
}

template <class ObjType, class RetType, class ValType1>
RetType simpleGetter (DynamicLibrary *lib, ObjType *obj,
	      std::string funcName, ValType1 parm1)
{ std::string errStr ;
  typedef RetType (*FuncSig)(ObjType *,ValType1) ;
  FuncSig getterFunc =
    lib->getFunc<FuncSig>(funcName,errStr) ;
  if (getterFunc != nullptr) {
    return (getterFunc(obj,parm1)) ;
  } else {
    std::cout
      << "  simpleGetter(1): lookup failed for \"" << funcName << "\": "
      << errStr << std::endl ;
  }
  return (RetType()) ;
}

template <class ObjType, class RetType, class ValType1,
	  class ValType2 = ValType1>
RetType simpleGetter (DynamicLibrary *lib, ObjType *obj,
	      std::string funcName, ValType1 parm1, ValType2 parm2)
{ std::string errStr ;
  typedef RetType (*FuncSig)(ObjType *,ValType1,ValType2) ;
  FuncSig getterFunc =
    lib->getFunc<FuncSig>(funcName,errStr) ;
  if (getterFunc != nullptr) {
    return (getterFunc(obj,parm1,parm2)) ;
  } else {
    std::cout
      << "  simpleGetter(2): lookup failed for \"" << funcName << "\": "
      << errStr << std::endl ;
  }
  return (RetType()) ;
}
  
template <class ObjType>
void simpleSetter (DynamicLibrary *lib, ObjType *obj, std::string funcName)
{ std::string errStr ;
  typedef void (*FuncSig)(ObjType *) ;
  FuncSig setterFunc =
    lib->getFunc<FuncSig>(funcName,errStr) ;
  if (setterFunc != nullptr) {
    setterFunc(obj) ;
  } else {
    std::cout
      << "  simpleSetter(0): lookup failed for \"" << funcName << "\": "
      << errStr << std::endl ;
  }
}

template <class ObjType, class ValType1>
void simpleSetter (DynamicLibrary *lib, ObjType *obj,
	      std::string funcName, ValType1 val1)
{ std::string errStr ;
  typedef void (*FuncSig)(ObjType *,ValType1) ;
  FuncSig setterFunc =
    lib->getFunc<FuncSig>(funcName,errStr) ;
  if (setterFunc != nullptr) {
    setterFunc(obj,val1) ;
  } else {
    std::cout
      << "  simpleSetter(1): lookup failed for \"" << funcName << "\": "
      << errStr << std::endl ;
  }
}

template <class ObjType, class ValType1, class ValType2 = ValType1>
void simpleSetter (DynamicLibrary *lib, ObjType *obj,
	      std::string funcName, ValType1 val1, ValType2 val2)
{ std::string errStr ;
  typedef void (*FuncSig)(ObjType *,ValType1,ValType2) ;
  FuncSig setterFunc =
    lib->getFunc<FuncSig>(funcName,errStr) ;
  if (setterFunc != nullptr) {
    setterFunc(obj,val1,val2) ;
  } else {
    std::cout
      << "  simpleSetter(2): lookup failed for \"" << funcName << "\": "
      << errStr << std::endl ;
  }
}

template <class ObjType,
	  class ValType1, class ValType2 = ValType1, class ValType3 = ValType1 >
void simpleSetter (DynamicLibrary *lib, ObjType *obj,
	    std::string funcName, ValType1 val1, ValType2 val2, ValType3 val3)
{ std::string errStr ;
  typedef void (*FuncSig)(ObjType *,ValType1,ValType2,ValType3) ;
  FuncSig setterFunc =
    lib->getFunc<FuncSig>(funcName,errStr) ;
  if (setterFunc != nullptr) {
    setterFunc(obj,val1,val2,val3) ;
  } else {
    std::cout
      << "  simpleSetter(3): lookup failed for \"" << funcName << "\": "
      << errStr << std::endl ;
  }
}

/*
  A pair of methods that assume the following model: There's a map of known
  members, keyed on a string that names some member of an object. Each entry
  in the map holds the names of the C get and set methods for the member. A
  get method takes no parameters other than the object and returns a value. A
  set method takes the object and a data value and returns void.

  `Member' is conceptual. The get and set methods can do pretty much anything
  as long as they match the signatures described above.

  In support, a typedef for the map and a macro to make it easy to create
  the map entries.
*/

#define mapEntry3(zzKey,zzGetName,zzSetName) \
{ std::string(zzKey), {std::string(zzGetName),std::string(zzSetName)} }

struct KMMFuncs { std::string getter_ ; std::string setter_ ; } ;
typedef std::map<std::string, KMMFuncs> KnownMemberMap ;


template <class ObjType, class ValType>
void setMember (KnownMemberMap membMap, std::string memb,
		DynamicLibrary *lib, ObjType *obj, ValType val)
{
  KnownMemberMap::const_iterator pxmIter = membMap.find(memb) ;
  if (pxmIter == membMap.end()) {
    std::cout
      << " Member " << memb << " is not recognised."
      << std::endl ;
    return ;
  }
  std::string setFuncName = pxmIter->second.setter_ ;
  simpleSetter<ObjType,ValType>(lib,obj,setFuncName,val) ;
}

template <class ObjType, class RetType>
RetType getMember (KnownMemberMap membMap, std::string memb,
		   DynamicLibrary *lib, ObjType *obj)
{
  KnownMemberMap::const_iterator pxmIter = membMap.find(memb) ;
  if (pxmIter == membMap.end()) {
    std::cout
      << " Member " << memb << " is not recognised."
      << std::endl ;
    return (RetType()) ;
  }
  std::string getFuncName = pxmIter->second.getter_ ;
  RetType retval = simpleGetter<ObjType,RetType>(lib,obj,getFuncName) ;
  return (retval) ;
}

