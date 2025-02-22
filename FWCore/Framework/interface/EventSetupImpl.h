#ifndef FWCore_Framework_EventSetupImpl_h
#define FWCore_Framework_EventSetupImpl_h
// -*- C++ -*-
//
// Package:     Framework
// Class:      EventSetupImpl
//
/**\class edm::EventSetupImpl

 Description: Container for all Records dealing with non-RunState info

 Usage:
    <usage>

*/
//
// Author:      Chris Jones
// Created:     Thu Mar 24 13:50:04 EST 2005
//

// system include files
#include <map>
#include <optional>
#include <vector>
#include "tbb/task_arena.h"

// user include files
#include "FWCore/Framework/interface/EventSetupRecordKey.h"
#include "FWCore/Framework/interface/EventSetupRecord.h"
#include "FWCore/Utilities/interface/ESIndices.h"

// forward declarations

class testEventsetup;

namespace edm {
  class ESInputTag;
  class ProcessBlockTransitionInfo;
  class Schedule;
  class ServiceToken;
  class WaitingTaskHolder;

  namespace eventsetup {
    class EventSetupProvider;
    class EventSetupRecordImpl;
    class EventSetupRecordProvider;
  }  // namespace eventsetup

  class EventSetupImpl {
  public:
    ~EventSetupImpl();
    EventSetupImpl() = delete;
    EventSetupImpl(EventSetupImpl const&) = delete;
    EventSetupImpl& operator=(EventSetupImpl const&) = delete;

    // ---------- const member functions ---------------------
    eventsetup::EventSetupRecordImpl const* findImpl(const eventsetup::EventSetupRecordKey&) const;
    eventsetup::EventSetupRecordImpl const* findImpl(ESRecordIndex) const;

    std::optional<eventsetup::EventSetupRecordGeneric> find(const eventsetup::EventSetupRecordKey&,
                                                            unsigned int iTransitionID,
                                                            ESProxyIndex const* getTokenIndices,
                                                            ESParentContext const& iParent) const;

    ///clears the oToFill vector and then fills it with the keys for all available records
    void fillAvailableRecordKeys(std::vector<eventsetup::EventSetupRecordKey>& oToFill) const;

    ///returns true if the Record is provided by a Source or a Producer
    /// a value of true does not mean this EventSetup object holds such a record
    bool recordIsProvidedByAModule(eventsetup::EventSetupRecordKey const&) const;

    bool validRecord(eventsetup::EventSetupRecordKey const& iKey) const;

    tbb::task_arena* taskArena CMS_THREAD_SAFE() const { return taskArena_; }
    ///Only EventSetupProvider allowed to create an EventSetupImpl
    friend class eventsetup::EventSetupProvider;
    friend class eventsetup::EventSetupRecordProvider;
    friend class ::testEventsetup;
    friend class ::testEventsetupRecord;
    friend class ProcessBlockTransitionInfo;

  protected:
    void addRecordImpl(const eventsetup::EventSetupRecordImpl& iRecord);

  private:
    explicit EventSetupImpl(tbb::task_arena*);

    void insertRecordImpl(const eventsetup::EventSetupRecordKey&, const eventsetup::EventSetupRecordImpl*);

    void setKeyIters(std::vector<eventsetup::EventSetupRecordKey>::const_iterator const& keysBegin,
                     std::vector<eventsetup::EventSetupRecordKey>::const_iterator const& keysEnd);

    // ---------- member data --------------------------------

    std::vector<eventsetup::EventSetupRecordKey>::const_iterator keysBegin_;
    std::vector<eventsetup::EventSetupRecordKey>::const_iterator keysEnd_;
    std::vector<eventsetup::EventSetupRecordImpl const*> recordImpls_;
    tbb::task_arena* taskArena_ = nullptr;
  };
}  // namespace edm
#endif
