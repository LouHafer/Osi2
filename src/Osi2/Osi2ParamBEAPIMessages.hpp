/*
  Copyright 2019 Lou Hafer, Matt Saltzman
  This code is licensed under the terms of the Eclipse Public License (EPL)
*/
/*! \file Osi2ParamBEAPIMessages.hpp
    \brief Type definitions for Osi2::ParamMgmtAPI messages.
*/

#ifndef Osi2ParamBEAPIMessages_HPP
# define Osi2ParamBEAPIMessages_HPP

#include "CoinMessageHandler.hpp"
#include "Osi2Config.h"

namespace Osi2 {

/*! \brief Enum to specify individual Osi2::ParamBEAPI messages

  The only purpose here is to define symbolic names for the messages. This
  makes it easy to add and remove messages.
*/

enum ParamBEAPIMsg {
    PMBEAPI_INIT,
    PMBEAPI_DUMMY_END
};

/*! \brief Osi2::ParamBEAPI messages

  This class holds the array of messages once they are loaded and provides the
  methods to use them. Its only reason for existence is the custom constructor,
  which loads the messages in the order defined by ParamMgmtAPIMsg. All else is
  inherited from CoinMessages.

  See CoinMessageHandler for usage.
*/
class OSI2LIB_EXPORT ParamBEAPIMessages : public CoinMessages {

public:

    /// Custom constructor
    ParamBEAPIMessages(Language language = us_en) ;

} ;

}  // end namespace Osi2

#endif

