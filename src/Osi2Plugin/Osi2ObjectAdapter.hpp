#ifndef OSI2OBJECT_ADAPTER_H
#define OSI2OBJECT_ADAPTER_H

#include "Osi2Plugin.hpp"

namespace Osi2 {

/*! \brief This interface is used to adapt C plugin objects to C++ plugin
	   objects.

  It must be passed to the PluginManager::createObject() function.
*/
struct IObjectAdapter {
    virtual ~IObjectAdapter() {}
    virtual void *adapt(void *object, DestroyFunc df) = 0 ;
} ;

/*! \brief This template should be used if the object model implements the
	   dual C/C++ object design pattern.

  Otherwise you need to provide your own object adapter class that implements
  IObjectAdapter

  /todo Seems to me this needs further implementation to work. In particular,
  what's the implementation of T((U *)object,df)? Seems like somewhere we need
  to provide a class definition for T. It needs to provide a constructor that
  accepts the object, a destructor that invokes DestroyFunc, and some set of
  methods that know how to invoke object (assuming object is a C function of
  some sort). Since we haven't grappled with a C plugin just yet, it's all
  speculation.  -- lh, 190711 --
*/
template <typename T, typename U>
struct ObjectAdapter : public IObjectAdapter {
    virtual void *adapt (void *object, DestroyFunc df) {
        return new T((U *)object, df);
    }
} ;


/*! \brief Dummy adapter

  For C++ plugins, nothing needs to happen to the returned object, and the
  adapted object isn't used. Still, there needs to be an implementation
  of a method that takes a void *object and a DestroyFunc. This is it.
*/

struct DummyAdapter : public IObjectAdapter {
  virtual void *adapt(void *object, DestroyFunc df)
  { return (static_cast<void*>(0)) ; }
} ;

}  // end namespace Osi2

#endif // OSI2OBJECT_ADAPTER_H
