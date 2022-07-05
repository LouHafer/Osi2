/*
  Copyright 2011 Lou Hafer, Matt Saltzman
  This code is licensed under the terms of the Eclipse Public License (EPL)
*/

#include <cstring>

#include "Osi2Config.h"

#include "Osi2CtrlAPIMessages.hpp"

namespace Osi2 {

/// Convenience structure for message definition.

typedef struct {
    /// Internal ID; used to specify a message
    CtrlAPIMsg intID_ ;

    /*! \brief External ID; printed in the message prefix

      Also tied to severity; see CoinMessageHandler
    */
    int extID_ ;

    /*! Message log level

      Message will print if log level is less than or equal to the current
      log level.
    */
    char logLvl_ ;

    /// Message text (may contain format codes).
    const char *text_ ;
} OneCtrlAPIMessage ;

/*
  Definitions of messages.  To add a new message:
    * Choose a unique ID. Add it to the enum in Osi2CtrlAPIMessages.hpp
    * Choose an external ID number from the appropriate range (information,
      warning, nonfatal or fatal error).
    * Define the message in the appropriate place in the list by adding the
      initialisation expression for a OneCtrlAPIMessage structure.

  Don't put anything after DUMMY_END. It's a marker used in the constructor.
*/
static OneCtrlAPIMessage us_english[] = {

    // Information: 0 -- 2999

    { CTRLAPI_INIT, 0, 7, "Control API constructor (%s)." },
    { CTRLAPI_LIBLDOK, 1, 7, "Plugin library \"%s\" (\"%s\") loaded." },
    { CTRLAPI_LIBCLOSEOK, 2, 7, "Plugin library \"%s\" (\"%s\") unloaded." },
    { CTRLAPI_CREATEOK, 3, 7, "API \"%s\"%?, library \"%s\"%? created." },
    { CTRLAPI_DESTROYOK, 4, 7, "API \"%s\"%?, library \"%s\"%? destroyed." },

    // Warning: 3000 -- 5999

    { CTRLAPI_LIBUNREG, 3000, 4, "Library \"%s\" is not registered." },
    { CTRLAPI_UNREG, 3001, 4,
        "%?PluginManager says \"%s\" already loaded but %?\"%s\" is not registered."
    },

    // Nonfatal Error: 6000 -- 8999

    {
        CTRLAPI_LIBLDFAIL, 6000, 2,
        "Load failed for plugin library \"%s\" (\"%s\")."
    },
    {
        CTRLAPI_LIBCLOSEFAIL, 6001, 2,
        "Load failed for plugin library \"%s\" (\"%s\")."
    },
    {
        CTRLAPI_CREATEFAIL, 6002, 2,
        "Create failed for API \"%s\"%?, library \"%s\"%?."
    },
    {
        CTRLAPI_DESTROYFAIL, 6003, 2,
        "Destroy failed for API \"%s\"%?, library \"%s\"%?."
    },
    { CTRLAPI_NOAPIIDENT, 6004, 2,
        "API object has no identity information!" },

    // Fatal Error: 9000 -- 9999

    { CTRLAPI_NOPLUGMGR, 9000, 1, "Cannot find plugin manager!" },

    { CTRLAPI_DUMMY_END, 9999, 0, "" }
} ;

/*
  Custom constructor. Really all we're doing here is loading each
  OneCtrlAPIMessage, in order according to internal ID, into a CoinOneMessage
  and installing it in the messages array.

  class_ is a gross hack. Arguably, we're class 2. See CoinMessageHandler.hpp.

  The string stored in source_ is limited to four chars (plus terminating
  null).
*/
CtrlAPIMessages::CtrlAPIMessages (Language language)
    : CoinMessages(sizeof(us_english) / sizeof(OneCtrlAPIMessage))
{
    language_ = language ;
    std::strcpy(source_, "OSI2") ;
    class_ = 2 ;

    OneCtrlAPIMessage *msg = &us_english[0] ;
    while (msg->intID_ != CTRLAPI_DUMMY_END) {
        CoinOneMessage tmpMsg(msg->extID_, msg->logLvl_, msg->text_) ;
        addMessage(msg->intID_, tmpMsg) ;
        msg++ ;
    }
    /*
      If we want to allow language-specific overrides of individual messages,
      code to do that should go here, prior to compacting the messages.
    */
    toCompact() ;
}

}  // end namespace Osi2

