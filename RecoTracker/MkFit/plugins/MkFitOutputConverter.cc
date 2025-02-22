#include "FWCore/Framework/interface/global/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/do_nothing_deleter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/TrackReco/interface/SeedStopInfo.h"
#include "DataFormats/TrackingRecHit/interface/InvalidTrackingRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"

#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "TrackingTools/TrackFitters/interface/KFTrajectoryFitter.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TkClonerImpl.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"

#include "RecoTracker/MkFit/interface/MkFitEventOfHits.h"
#include "RecoTracker/MkFit/interface/MkFitClusterIndexToHit.h"
#include "RecoTracker/MkFit/interface/MkFitSeedWrapper.h"
#include "RecoTracker/MkFit/interface/MkFitOutputWrapper.h"
#include "RecoTracker/MkFit/interface/MkFitGeometry.h"
#include "RecoTracker/Record/interface/TrackerRecoGeometryRecord.h"

// mkFit indludes
#include "LayerNumberConverter.h"
#include "Track.h"
#include "mkFit/HitStructures.h"

namespace {
  template <typename T>
  bool isBarrel(T subdet) {
    return subdet == PixelSubdetector::PixelBarrel || subdet == StripSubdetector::TIB ||
           subdet == StripSubdetector::TOB;
  }

  template <typename T>
  bool isEndcap(T subdet) {
    return subdet == PixelSubdetector::PixelEndcap || subdet == StripSubdetector::TID ||
           subdet == StripSubdetector::TEC;
  }
}  // namespace

class MkFitOutputConverter : public edm::global::EDProducer<> {
public:
  explicit MkFitOutputConverter(edm::ParameterSet const& iConfig);
  ~MkFitOutputConverter() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const override;

  TrackCandidateCollection convertCandidates(const MkFitOutputWrapper& mkFitOutput,
                                             const mkfit::EventOfHits& eventOfHits,
                                             const MkFitClusterIndexToHit& pixelClusterIndexToHit,
                                             const MkFitClusterIndexToHit& stripClusterIndexToHit,
                                             const edm::View<TrajectorySeed>& seeds,
                                             const MagneticField& mf,
                                             const Propagator& propagatorAlong,
                                             const Propagator& propagatorOpposite,
                                             const TkClonerImpl& hitCloner,
                                             const std::vector<const DetLayer*>& detLayers,
                                             const mkfit::TrackVec& mkFitSeeds) const;

  std::pair<TrajectoryStateOnSurface, const GeomDet*> backwardFit(const FreeTrajectoryState& fts,
                                                                  const edm::OwnVector<TrackingRecHit>& hits,
                                                                  const Propagator& propagatorAlong,
                                                                  const Propagator& propagatorOpposite,
                                                                  const TkClonerImpl& hitCloner,
                                                                  bool lastHitWasInvalid,
                                                                  bool lastHitWasChanged) const;

  std::pair<TrajectoryStateOnSurface, const GeomDet*> convertInnermostState(const FreeTrajectoryState& fts,
                                                                            const edm::OwnVector<TrackingRecHit>& hits,
                                                                            const Propagator& propagatorAlong,
                                                                            const Propagator& propagatorOpposite) const;

  const edm::EDGetTokenT<MkFitEventOfHits> eventOfHitsToken_;
  const edm::EDGetTokenT<MkFitClusterIndexToHit> pixelClusterIndexToHitToken_;
  const edm::EDGetTokenT<MkFitClusterIndexToHit> stripClusterIndexToHitToken_;
  const edm::EDGetTokenT<MkFitSeedWrapper> mkfitSeedToken_;
  const edm::EDGetTokenT<MkFitOutputWrapper> tracksToken_;
  const edm::EDGetTokenT<edm::View<TrajectorySeed>> seedToken_;
  const edm::ESGetToken<Propagator, TrackingComponentsRecord> propagatorAlongToken_;
  const edm::ESGetToken<Propagator, TrackingComponentsRecord> propagatorOppositeToken_;
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> mfToken_;
  const edm::ESGetToken<TransientTrackingRecHitBuilder, TransientRecHitRecord> ttrhBuilderToken_;
  const edm::ESGetToken<MkFitGeometry, TrackerRecoGeometryRecord> mkFitGeomToken_;
  const edm::EDPutTokenT<TrackCandidateCollection> putTrackCandidateToken_;
  const edm::EDPutTokenT<std::vector<SeedStopInfo>> putSeedStopInfoToken_;
  const std::string ttrhBuilderName_;
  const std::string propagatorAlongName_;
  const std::string propagatorOppositeName_;
};

MkFitOutputConverter::MkFitOutputConverter(edm::ParameterSet const& iConfig)
    : eventOfHitsToken_{consumes<MkFitEventOfHits>(iConfig.getParameter<edm::InputTag>("mkFitEventOfHits"))},
      pixelClusterIndexToHitToken_{consumes(iConfig.getParameter<edm::InputTag>("mkFitPixelHits"))},
      stripClusterIndexToHitToken_{consumes(iConfig.getParameter<edm::InputTag>("mkFitStripHits"))},
      mkfitSeedToken_{consumes<MkFitSeedWrapper>(iConfig.getParameter<edm::InputTag>("mkFitSeeds"))},
      tracksToken_{consumes<MkFitOutputWrapper>(iConfig.getParameter<edm::InputTag>("tracks"))},
      seedToken_{consumes<edm::View<TrajectorySeed>>(iConfig.getParameter<edm::InputTag>("seeds"))},
      propagatorAlongToken_{
          esConsumes<Propagator, TrackingComponentsRecord>(iConfig.getParameter<edm::ESInputTag>("propagatorAlong"))},
      propagatorOppositeToken_{esConsumes<Propagator, TrackingComponentsRecord>(
          iConfig.getParameter<edm::ESInputTag>("propagatorOpposite"))},
      mfToken_{esConsumes<MagneticField, IdealMagneticFieldRecord>()},
      ttrhBuilderToken_{esConsumes<TransientTrackingRecHitBuilder, TransientRecHitRecord>(
          iConfig.getParameter<edm::ESInputTag>("ttrhBuilder"))},
      mkFitGeomToken_{esConsumes<MkFitGeometry, TrackerRecoGeometryRecord>()},
      putTrackCandidateToken_{produces<TrackCandidateCollection>()},
      putSeedStopInfoToken_{produces<std::vector<SeedStopInfo>>()} {}

void MkFitOutputConverter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add("mkFitEventOfHits", edm::InputTag{"mkFitEventOfHits"});
  desc.add("mkFitPixelHits", edm::InputTag{"mkFitSiPixelHits"});
  desc.add("mkFitStripHits", edm::InputTag{"mkFitSiStripHits"});
  desc.add("mkFitSeeds", edm::InputTag{"mkFitSeedConverter"});
  desc.add("tracks", edm::InputTag{"mkFitProducer"});
  desc.add("seeds", edm::InputTag{"initialStepSeeds"});
  desc.add("ttrhBuilder", edm::ESInputTag{"", "WithTrackAngle"});
  desc.add("propagatorAlong", edm::ESInputTag{"", "PropagatorWithMaterial"});
  desc.add("propagatorOpposite", edm::ESInputTag{"", "PropagatorWithMaterialOpposite"});

  descriptions.addWithDefaultLabel(desc);
}

void MkFitOutputConverter::produce(edm::StreamID iID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  const auto& seeds = iEvent.get(seedToken_);
  const auto& mkfitSeeds = iEvent.get(mkfitSeedToken_);

  const auto& ttrhBuilder = iSetup.getData(ttrhBuilderToken_);
  const auto* tkBuilder = dynamic_cast<TkTransientTrackingRecHitBuilder const*>(&ttrhBuilder);
  if (!tkBuilder) {
    throw cms::Exception("LogicError") << "TTRHBuilder must be of type TkTransientTrackingRecHitBuilder";
  }
  const auto& mkFitGeom = iSetup.getData(mkFitGeomToken_);

  // Convert mkfit presentation back to CMSSW
  iEvent.emplace(putTrackCandidateToken_,
                 convertCandidates(iEvent.get(tracksToken_),
                                   iEvent.get(eventOfHitsToken_).get(),
                                   iEvent.get(pixelClusterIndexToHitToken_),
                                   iEvent.get(stripClusterIndexToHitToken_),
                                   seeds,
                                   iSetup.getData(mfToken_),
                                   iSetup.getData(propagatorAlongToken_),
                                   iSetup.getData(propagatorOppositeToken_),
                                   tkBuilder->cloner(),
                                   mkFitGeom.detLayers(),
                                   mkfitSeeds.seeds()));

  // TODO: SeedStopInfo is currently unfilled
  iEvent.emplace(putSeedStopInfoToken_, seeds.size());
}

TrackCandidateCollection MkFitOutputConverter::convertCandidates(const MkFitOutputWrapper& mkFitOutput,
                                                                 const mkfit::EventOfHits& eventOfHits,
                                                                 const MkFitClusterIndexToHit& pixelClusterIndexToHit,
                                                                 const MkFitClusterIndexToHit& stripClusterIndexToHit,
                                                                 const edm::View<TrajectorySeed>& seeds,
                                                                 const MagneticField& mf,
                                                                 const Propagator& propagatorAlong,
                                                                 const Propagator& propagatorOpposite,
                                                                 const TkClonerImpl& hitCloner,
                                                                 const std::vector<const DetLayer*>& detLayers,
                                                                 const mkfit::TrackVec& mkFitSeeds) const {
  TrackCandidateCollection output;
  const auto& candidates = mkFitOutput.tracks();
  output.reserve(candidates.size());

  LogTrace("MkFitOutputConverter") << "Number of candidates " << candidates.size();

  int candIndex = -1;
  for (const auto& cand : candidates) {
    ++candIndex;
    LogTrace("MkFitOutputConverter") << "Candidate " << candIndex << " pT " << cand.pT() << " eta " << cand.momEta()
                                     << " phi " << cand.momPhi() << " chi2 " << cand.chi2();

    // hits
    edm::OwnVector<TrackingRecHit> recHits;
    // nTotalHits() gives sum of valid hits (nFoundHits()) and invalid/missing hits.
    const int nhits = cand.nTotalHits();
    bool lastHitInvalid = false;
    for (int i = 0; i < nhits; ++i) {
      const auto& hitOnTrack = cand.getHitOnTrack(i);
      LogTrace("MkFitOutputConverter") << " hit on layer " << hitOnTrack.layer << " index " << hitOnTrack.index;
      if (hitOnTrack.index < 0) {
        // See index-desc.txt file in mkFit for description of negative values
        //
        // In order to use the regular InvalidTrackingRecHit I'd need
        // a GeomDet (and "unfortunately" that is needed in
        // TrackProducer).
        //
        // I guess we could take the track state and propagate it to
        // each layer to find the actual module the track crosses, and
        // check whether it is active or not to be able to mark
        // inactive hits
        const auto* detLayer = detLayers.at(hitOnTrack.layer);
        if (detLayer == nullptr) {
          throw cms::Exception("LogicError") << "DetLayer for layer index " << hitOnTrack.layer << " is null!";
        }
        // In principle an InvalidTrackingRecHitNoDet could be
        // inserted here, but it seems that it is best to deal with
        // them in the TrackProducer.
        lastHitInvalid = true;
      } else {
        auto const isPixel = eventOfHits[hitOnTrack.layer].is_pix_lyr();
        auto const& hits = isPixel ? pixelClusterIndexToHit.hits() : stripClusterIndexToHit.hits();

        auto const& thit = static_cast<BaseTrackerRecHit const&>(*hits[hitOnTrack.index]);
        if (thit.firstClusterRef().isPixel() || thit.detUnit()->type().isEndcap()) {
          recHits.push_back(hits[hitOnTrack.index]->clone());
        } else {
          recHits.push_back(std::make_unique<SiStripRecHit1D>(
              thit.localPosition(),
              LocalError(thit.localPositionError().xx(), 0.f, std::numeric_limits<float>::max()),
              *thit.det(),
              thit.firstClusterRef()));
        }
        LogTrace("MkFitOutputConverter") << "  pos " << recHits.back().globalPosition().x() << " "
                                         << recHits.back().globalPosition().y() << " "
                                         << recHits.back().globalPosition().z() << " mag2 "
                                         << recHits.back().globalPosition().mag2() << " detid "
                                         << recHits.back().geographicalId().rawId() << " cluster " << hitOnTrack.index;
        lastHitInvalid = false;
      }
    }

    const auto lastHitId = recHits.back().geographicalId();

    // MkFit hits are *not* in the order of propagation, sort by 3D radius for now (as we don't have loopers)
    // TODO: Improve the sorting (extract keys? maybe even bubble sort would work well as the hits are almost in the correct order)
    recHits.sort([](const auto& a, const auto& b) {
      const auto asub = a.geographicalId().subdetId();
      const auto bsub = b.geographicalId().subdetId();
      if (asub != bsub) {
        // Subdetector order (BPix, FPix, TIB, TID, TOB, TEC) corresponds also the navigation
        return asub < bsub;
      }

      const auto& apos = a.globalPosition();
      const auto& bpos = b.globalPosition();

      if (isBarrel(asub)) {
        return apos.perp2() < bpos.perp2();
      }
      return std::abs(apos.z()) < std::abs(bpos.z());
    });

    const bool lastHitChanged = (recHits.back().geographicalId() != lastHitId);  // TODO: make use of the bools

    // seed
    const auto seedIndex = cand.label();
    LogTrace("MkFitOutputConverter") << " from seed " << seedIndex << " seed hits";

    // state
    auto state = cand.state();  // copy because have to modify
    state.convertFromCCSToGlbCurvilinear();
    const auto& param = state.parameters;
    const auto& err = state.errors;
    AlgebraicSymMatrix55 cov;
    for (int i = 0; i < 5; ++i) {
      for (int j = i; j < 5; ++j) {
        cov[i][j] = err.At(i, j);
      }
    }

    auto fts = FreeTrajectoryState(
        GlobalTrajectoryParameters(
            GlobalPoint(param[0], param[1], param[2]), GlobalVector(param[3], param[4], param[5]), state.charge, &mf),
        CurvilinearTrajectoryError(cov));
    if (!fts.curvilinearError().posDef()) {
      edm::LogWarning("MkFitOutputConverter") << "Curvilinear error not pos-def\n"
                                              << fts.curvilinearError().matrix() << "\ncandidate ignored";
      continue;
    }

    auto tsosDet =
        mkFitOutput.propagatedToFirstLayer()
            ? convertInnermostState(fts, recHits, propagatorAlong, propagatorOpposite)
            : backwardFit(fts, recHits, propagatorAlong, propagatorOpposite, hitCloner, lastHitInvalid, lastHitChanged);
    if (!tsosDet.first.isValid()) {
      edm::LogWarning("MkFitOutputConverter")
          << "Backward fit of candidate " << candIndex << " failed, ignoring the candidate";
      continue;
    }

    // convert to persistent, from CkfTrackCandidateMakerBase
    auto pstate = trajectoryStateTransform::persistentState(tsosDet.first, tsosDet.second->geographicalId().rawId());

    output.emplace_back(
        recHits,
        seeds.at(seedIndex),
        pstate,
        seeds.refAt(seedIndex),
        0,                                               // mkFit does not produce loopers, so set nLoops=0
        static_cast<uint8_t>(StopReason::UNINITIALIZED)  // TODO: ignore details of stopping reason as well for now
    );
  }
  return output;
}

std::pair<TrajectoryStateOnSurface, const GeomDet*> MkFitOutputConverter::backwardFit(
    const FreeTrajectoryState& fts,
    const edm::OwnVector<TrackingRecHit>& hits,
    const Propagator& propagatorAlong,
    const Propagator& propagatorOpposite,
    const TkClonerImpl& hitCloner,
    bool lastHitWasInvalid,
    bool lastHitWasChanged) const {
  // First filter valid hits as in TransientInitialStateEstimator
  TransientTrackingRecHit::ConstRecHitContainer firstHits;

  for (int i = hits.size() - 1; i >= 0; --i) {
    if (hits[i].det()) {
      // TransientTrackingRecHit::ConstRecHitContainer has shared_ptr,
      // and it is passed to backFitter below so it is really needed
      // to keep the interface. Since we keep the ownership in hits,
      // let's disable the deleter.
      firstHits.emplace_back(&(hits[i]), edm::do_nothing_deleter{});
    }
  }

  // Then propagate along to the surface of the last hit to get a TSOS
  const auto& lastHitSurface = firstHits.front()->det()->surface();

  const Propagator* tryFirst = &propagatorAlong;
  const Propagator* trySecond = &propagatorOpposite;
  if (lastHitWasInvalid || lastHitWasChanged) {
    LogTrace("MkFitOutputConverter") << "Propagating first opposite, then along, because lastHitWasInvalid? "
                                     << lastHitWasInvalid << " or lastHitWasChanged? " << lastHitWasChanged;
    std::swap(tryFirst, trySecond);
  } else {
    const auto lastHitSubdet = firstHits.front()->geographicalId().subdetId();
    const auto& surfacePos = lastHitSurface.position();
    const auto& lastHitPos = firstHits.front()->globalPosition();
    bool doSwitch = false;
    if (isBarrel(lastHitSubdet)) {
      doSwitch = (surfacePos.perp2() < lastHitPos.perp2());
    } else {
      doSwitch = (surfacePos.z() < lastHitPos.z());
    }
    if (doSwitch) {
      LogTrace("MkFitOutputConverter")
          << "Propagating first opposite, then along, because surface is inner than the hit; surface perp2 "
          << surfacePos.perp() << " hit " << lastHitPos.perp2() << " surface z " << surfacePos.z() << " hit "
          << lastHitPos.z();

      std::swap(tryFirst, trySecond);
    }
  }

  auto tsosDouble = tryFirst->propagateWithPath(fts, lastHitSurface);
  if (!tsosDouble.first.isValid()) {
    LogDebug("MkFitOutputConverter") << "Propagating to startingState failed, trying in another direction next";
    tsosDouble = trySecond->propagateWithPath(fts, lastHitSurface);
  }
  auto& startingState = tsosDouble.first;

  if (!startingState.isValid()) {
    edm::LogWarning("MkFitOutputConverter")
        << "startingState is not valid, FTS was\n"
        << fts << " last hit surface surface:"
        << "\n position " << lastHitSurface.position() << "\n phiSpan " << lastHitSurface.phiSpan().first << ","
        << lastHitSurface.phiSpan().first << "\n rSpan " << lastHitSurface.rSpan().first << ","
        << lastHitSurface.rSpan().first << "\n zSpan " << lastHitSurface.zSpan().first << ","
        << lastHitSurface.zSpan().first;
    return std::pair<TrajectoryStateOnSurface, const GeomDet*>();
  }

  // Then return back to the logic from TransientInitialStateEstimator
  startingState.rescaleError(100.);

  // avoid cloning
  KFUpdator const aKFUpdator;
  Chi2MeasurementEstimator const aChi2MeasurementEstimator(100., 3);
  KFTrajectoryFitter backFitter(
      &propagatorAlong, &aKFUpdator, &aChi2MeasurementEstimator, firstHits.size(), nullptr, &hitCloner);

  // assume for now that the propagation in mkfit always alongMomentum
  PropagationDirection backFitDirection = oppositeToMomentum;

  // only direction matters in this context
  TrajectorySeed fakeSeed(PTrajectoryStateOnDet(), edm::OwnVector<TrackingRecHit>(), backFitDirection);

  // ignore loopers for now
  Trajectory fitres = backFitter.fitOne(fakeSeed, firstHits, startingState, TrajectoryFitter::standard);

  LogDebug("MkFitOutputConverter") << "using a backward fit of :" << firstHits.size() << " hits, starting from:\n"
                                   << startingState << " to get the estimate of the initial state of the track.";

  if (!fitres.isValid()) {
    edm::LogWarning("MkFitOutputConverter") << "FitTester: first hits fit failed";
    return std::pair<TrajectoryStateOnSurface, const GeomDet*>();
  }

  TrajectoryMeasurement const& firstMeas = fitres.lastMeasurement();

  // magnetic field can be different!
  TrajectoryStateOnSurface firstState(firstMeas.updatedState().localParameters(),
                                      firstMeas.updatedState().localError(),
                                      firstMeas.updatedState().surface(),
                                      propagatorAlong.magneticField());

  firstState.rescaleError(100.);

  LogDebug("MkFitOutputConverter") << "the initial state is found to be:\n:" << firstState
                                   << "\n it's field pointer is: " << firstState.magneticField()
                                   << "\n the pointer from the state of the back fit was: "
                                   << firstMeas.updatedState().magneticField();

  return std::make_pair(firstState, firstMeas.recHit()->det());
}

std::pair<TrajectoryStateOnSurface, const GeomDet*> MkFitOutputConverter::convertInnermostState(
    const FreeTrajectoryState& fts,
    const edm::OwnVector<TrackingRecHit>& hits,
    const Propagator& propagatorAlong,
    const Propagator& propagatorOpposite) const {
  auto det = hits[0].det();
  if (det == nullptr) {
    throw cms::Exception("LogicError") << "Got nullptr from the first hit det()";
  }

  const auto& firstHitSurface = det->surface();

  auto tsosDouble = propagatorAlong.propagateWithPath(fts, firstHitSurface);
  if (!tsosDouble.first.isValid()) {
    LogDebug("MkFitOutputConverter") << "Propagating to startingState along momentum failed, trying opposite next";
    tsosDouble = propagatorOpposite.propagateWithPath(fts, firstHitSurface);
  }

  return std::make_pair(tsosDouble.first, det);
}

DEFINE_FWK_MODULE(MkFitOutputConverter);
