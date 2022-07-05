/*
  Copyright 2011 Lou Hafer, Matt Saltzman
  This code is licensed under the terms of the Eclipse Public License (EPL)

  $Id$

  This file contains the implementation of the methods for the Osi2 ParamBE
  API.
*/

#include "Osi2Config.h"
#include "Osi2ParamBEAPI_Imp.hpp"
#include "Osi2ControlAPI.hpp"

namespace Osi2 {

template<class Client>
void ParamBEAPI_Imp<Client>::addParam (const char *paramID,
				        ParamEntry *entry)
{
  if (paramIDCnt_ >= paramIDLen_) {
    paramIDLen_ += 10 ;
    const char **tmp = new const char*[paramIDLen_] ;
    int bytes = paramIDLen_*sizeof(char *) ;
    for (int ndx = 0 ; ndx < paramIDCnt_ ; ndx++) {
      tmp[ndx] = paramIDs_[ndx] ;
    }
    delete[] paramIDs_ ;
    paramIDs_ = tmp ;
  }
  paramIDs_[paramIDCnt_++] = paramID ;
  entryForParam_[paramID] = entry ;

  return ;
}


/*
  Ask ParamEntry to invoke the appropriate get function. We already know the
  client's class. ParamEntry knows the parameter's type. 
*/
template<class Client>
bool ParamBEAPI_Imp<Client>::get (const char *paramID, void *&blob)
{
  ParamEntry *paramEntry = entryForParam_[paramID] ;
  paramEntry->get(client_,blob) ;
  return (true) ;
}


/*
  And a matching set method.
*/
template<class Client>
bool ParamBEAPI_Imp<Client>::set (const char *paramID, const void *&blob)
{
  ParamEntry *paramEntry = entryForParam_[paramID] ;
  paramEntry->set(client_,blob) ;
  return (true) ;
}

/*
  Boilerplate: Constructors, destructors, & such like
*/

/*
  Constructor. There's no default here; we need a client.
*/
template<class Client>
ParamBEAPI_Imp<Client>::ParamBEAPI_Imp (Client *client)
  : client_(client),
    paramIDCnt_(0),
    paramIDLen_(0),
    paramIDs_(nullptr),
    logLvl_(7)
{
  msgHandler_ = new CoinMessageHandler() ;
  msgs_ = ParamBEAPIMessages() ;
  dfltHandler_ = true ;
  msgHandler_->setLogLevel(logLvl_) ;
  msgHandler_->message(PMBEAPI_INIT, msgs_) << "default" << CoinMessageEol ;
}

/*
  Copy constructor
*/
template<class Client>
ParamBEAPI_Imp<Client>::ParamBEAPI_Imp (const ParamBEAPI_Imp &rhs)
  : client_(rhs.client_),
    apiMgr_(rhs.apiMgr_),
    dfltHandler_(rhs.dfltHandler_),
    logLvl_(rhs.logLvl_)
{
/*
  Similarly for the param ID vector.
*/
  paramIDCnt_ = rhs.paramIDCnt_ ;
  paramIDLen_ = rhs.paramIDLen_ ;
  paramIDs_ = new const char*[paramIDLen_] ;
  for (int ndx = 0 ; ndx < paramIDCnt_ ; ndx++) {
    paramIDs_[ndx] = rhs.paramIDs_[ndx] ;
  }
  for (typename EntryForParamMap::const_iterator iter =
  		rhs.entryForParam_.begin() ;
       iter != rhs.entryForParam_.end() ;
       iter++) {
    ParamEntry *entry = iter->second->clone() ;
    entryForParam_[iter->first] = entry ;
  }
/*
  If this is our handler, make an independent copy. If it's the client's
  handler, we can't make an independent copy because the client won't know
  about it and won't delete it.
*/
  if (dfltHandler_) {
    msgHandler_ = new CoinMessageHandler(*rhs.msgHandler_) ;
  } else {
    msgHandler_ = rhs.msgHandler_ ;
  }
  msgs_ = rhs.msgs_ ;
  msgHandler_->setLogLevel(logLvl_) ;
  msgHandler_->message(PMBEAPI_INIT, msgs_) << "copy" << CoinMessageEol ;
}

/*
  Assignment
*/
template<class Client>
ParamBEAPI_Imp<Client> &ParamBEAPI_Imp<Client>::operator=
    (const ParamBEAPI_Imp &rhs)
{
/*
  Self-assignment requires no work.
*/
  if (this == &rhs) return (*this) ;
/*
  Don't lose the client, or the APIs.
*/
  client_ = rhs.client_ ;
  apiMgr_ = rhs.apiMgr_ ;
/*
  Similarly for the param ID vector. The entryForParm_ map is a different
  story. It holds pointers to ParamEntry objects that we need to replicate.
*/
  paramIDCnt_ = rhs.paramIDCnt_ ;
  paramIDLen_ = rhs.paramIDLen_ ;
  paramIDs_ = new const char*[paramIDLen_] ;
  for (int ndx = 0 ; ndx < paramIDCnt_ ; ndx++) {
    paramIDs_[ndx] = rhs.paramIDs_[ndx] ;
  }
  for (typename EntryForParamMap::const_iterator iter =
  		    rhs.entryForParam_.begin() ;
       iter != rhs.entryForParam_.end() ;
       iter++) {
    ParamEntry *entry = iter->second->clone() ;
    entryForParam_[iter->first] = entry ;
  }
/*
  If it's our handler, we need to delete the old and replace with the new.
  If it's the user's handler, it's the user's problem. We just assign the
  pointer.
*/
  if (dfltHandler_) {
    delete msgHandler_ ;
    msgHandler_ = nullptr ;
  }
  dfltHandler_ = rhs.dfltHandler_ ;
  if (dfltHandler_) {
    msgHandler_ = new CoinMessageHandler(*rhs.msgHandler_) ;
  } else {
    msgHandler_ = rhs.msgHandler_ ;
  }
  msgs_ = rhs.msgs_ ;
  msgHandler_->setLogLevel(logLvl_) ;

  return (*this) ;
}


/*
  Destructor
*/
template<class Client>
ParamBEAPI_Imp<Client>::~ParamBEAPI_Imp ()
{
/*
  If we're holding a param ID vector, delete it. We're not responsible for the
  character strings it points to, just the vector.
*/
  if (paramIDs_) { delete[] paramIDs_ ; }
/*
  We are responsible for the ParamEntry objects held in entryForParam_
*/
  typedef typename std::map<std::string,ParamEntry *>::iterator EFPIter ;
  for (EFPIter iter = entryForParam_.begin() ;
       iter != entryForParam_.end() ;
       iter++) {
    if (iter->second) { delete iter->second ; }
  }

/*
  If this is our handler, delete it. Otherwise it's the client's
  responsibility.
*/
  if (dfltHandler_) {
    delete msgHandler_ ;
    msgHandler_ = nullptr ;
  }
}

/*
  Virtual constructor
*/
template<class Client>
ParamBEAPI *ParamBEAPI_Imp<Client>::create (Client *client)
{
  ParamBEAPI *api = new ParamBEAPI_Imp(client) ;
  return (api) ;
}

/*
  Clone
*/
template<class Client>
ParamBEAPI *ParamBEAPI_Imp<Client>::clone ()
{
    ParamBEAPI *api = new ParamBEAPI_Imp(*this) ;
    return (api) ;
}


}    // namespace Osi2
