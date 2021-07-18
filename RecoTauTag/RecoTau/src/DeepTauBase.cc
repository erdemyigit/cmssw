/*
 * \class DeepTauBase
 *
 * Implementation of the base class for tau identification using Deep NN.
 *
 * \author Konstantin Androsov, INFN Pisa
 * \author Maria Rosaria Di Domenico, University of Siena & INFN Pisa
 */

//TODO: port to offline RECO/AOD inputs to allow usage with offline AOD
//TODO: Take into account that PFTaus can also be build with pat::PackedCandidates

#include "RecoTauTag/RecoTau/interface/DeepTauBase.h"

namespace deep_tau {

  TauWPThreshold::TauWPThreshold(const std::string& cut_str) {
    bool simple_value = false;
    try {
      size_t pos = 0;
      value_ = std::stod(cut_str, &pos);
      simple_value = (pos == cut_str.size());
    } catch (std::invalid_argument&) {
    } catch (std::out_of_range&) {
    }
    if (!simple_value) {
      static const std::string prefix =
          "[&](double *x, double *p) { const int decayMode = p[0];"
          "const double pt = p[1]; const double eta = p[2];";
      static const int n_params = 3;
      static const auto handler = [](int, Bool_t, const char*, const char*) -> void {};

      const std::string fn_str = prefix + cut_str + "}";
      auto old_handler = SetErrorHandler(handler);
      fn_ = std::make_unique<TF1>("fn_", fn_str.c_str(), 0, 1, n_params);
      SetErrorHandler(old_handler);
      if (!fn_->IsValid())
        throw cms::Exception("TauWPThreshold: invalid formula") << "Invalid WP cut formula = '" << cut_str << "'.";
    }
  }

  double TauWPThreshold::operator()(const reco::BaseTau& tau, bool isPFTau) const {
    if (!fn_)
      return value_;

    if (isPFTau)
      fn_->SetParameter(0, dynamic_cast<const reco::PFTau&>(tau).decayMode());
    else
      fn_->SetParameter(0, dynamic_cast<const pat::Tau&>(tau).decayMode());
    fn_->SetParameter(1, tau.pt());
    fn_->SetParameter(2, tau.eta());
    return fn_->Eval(0);
  }

  float DeepTauBase::Output::get_number(const tensorflow::Tensor& pred, size_t tau_index, size_t elem) const {
    return pred.matrix<float>()(tau_index, elem);
  }

  float DeepTauBase::Output::get_number(const std::vector<std::vector<float>>& pred,
                                        size_t tau_index,
                                        size_t elem) const {
    return pred.at(tau_index).at(elem);
  }

  DeepTauBase::DeepTauBase(const edm::ParameterSet& cfg,
                           const OutputCollection& outputCollection,
                           const DeepTauCache* cache)
      : tausToken_(consumes<TauCollection>(cfg.getParameter<edm::InputTag>("taus"))),
        pfcandToken_(consumes<CandidateCollection>(cfg.getParameter<edm::InputTag>("pfcands"))),
        vtxToken_(consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("vertices"))),
        is_online_(cfg.getParameter<bool>("is_online")),
        outputs_(outputCollection),
        cache_(cache) {
    for (const auto& output_desc : outputs_) {
      produces<TauDiscriminator>(output_desc.first);
      const auto& cut_list = cfg.getParameter<std::vector<std::string>>(output_desc.first + "WP");
      for (const std::string& cut_str : cut_list) {
        workingPoints_[output_desc.first].push_back(std::make_unique<Cutter>(cut_str));
      }
    }

    // prediscriminant operator
    // require the tau to pass the following prediscriminants
    const edm::ParameterSet& prediscriminantConfig = cfg.getParameter<edm::ParameterSet>("Prediscriminants");

    // determine boolean operator used on the prediscriminants
    std::string pdBoolOperator = prediscriminantConfig.getParameter<std::string>("BooleanOperator");
    // convert string to lowercase
    transform(pdBoolOperator.begin(), pdBoolOperator.end(), pdBoolOperator.begin(), ::tolower);

    if (pdBoolOperator == "and") {
      andPrediscriminants_ = 0x1;  //use chars instead of bools so we can do a bitwise trick later
    } else if (pdBoolOperator == "or") {
      andPrediscriminants_ = 0x0;
    } else {
      throw cms::Exception("TauDiscriminationProducerBase")
          << "PrediscriminantBooleanOperator defined incorrectly, options are: AND,OR";
    }

    // get the list of prediscriminants
    std::vector<std::string> prediscriminantsNames =
        prediscriminantConfig.getParameterNamesForType<edm::ParameterSet>();

    for (auto const& iDisc : prediscriminantsNames) {
      const edm::ParameterSet& iPredisc = prediscriminantConfig.getParameter<edm::ParameterSet>(iDisc);
      const edm::InputTag& label = iPredisc.getParameter<edm::InputTag>("Producer");
      double cut = iPredisc.getParameter<double>("cut");

      if (is_online_) {
        TauDiscInfo<reco::PFTauDiscriminator> thisDiscriminator;
        thisDiscriminator.label = label;
        thisDiscriminator.cut = cut;
        thisDiscriminator.disc_token = consumes<reco::PFTauDiscriminator>(label);
        recoPrediscriminants_.push_back(thisDiscriminator);
      } else {
        TauDiscInfo<pat::PATTauDiscriminator> thisDiscriminator;
        thisDiscriminator.label = label;
        thisDiscriminator.cut = cut;
        thisDiscriminator.disc_token = consumes<pat::PATTauDiscriminator>(label);
        patPrediscriminants_.push_back(thisDiscriminator);
      }
    }
  }

  void DeepTauBase::produce(edm::Event& event, const edm::EventSetup& es) {
    edm::Handle<TauCollection> taus;
    event.getByToken(tausToken_, taus);
    edm::ProductID tauProductID = taus.id();

    // load prediscriminators
    size_t nPrediscriminants =
        patPrediscriminants_.empty() ? recoPrediscriminants_.size() : patPrediscriminants_.size();
    for (size_t iDisc = 0; iDisc < nPrediscriminants; ++iDisc) {
      edm::ProductID discKeyId;
      if (is_online_) {
        recoPrediscriminants_[iDisc].fill(event);
        discKeyId = recoPrediscriminants_[iDisc].handle->keyProduct().id();
      } else {
        patPrediscriminants_[iDisc].fill(event);
        discKeyId = patPrediscriminants_[iDisc].handle->keyProduct().id();
      }

      // Check to make sure the product is correct for the discriminator.
      // If not, throw a more informative exception.
      if (tauProductID != discKeyId) {
        throw cms::Exception("MisconfiguredPrediscriminant")
            << "The tau collection has product ID: " << tauProductID
            << " but the pre-discriminator is keyed with product ID: " << discKeyId << std::endl;
      }
    }

    const tensorflow::Tensor& pred = getPredictions(event, taus);
    createOutputs(event, pred, taus);
  }

  void DeepTauBase::createOutputs(edm::Event& event, const tensorflow::Tensor& pred, edm::Handle<TauCollection> taus) {
    for (const auto& output_desc : outputs_) {
      const WPList* working_points = nullptr;
      if (workingPoints_.find(output_desc.first) != workingPoints_.end()) {
        working_points = &workingPoints_.at(output_desc.first);
      }
      auto result = output_desc.second.get_value(taus, pred, working_points, is_online_);
      event.put(std::move(result), output_desc.first);
    }
  }

  std::unique_ptr<DeepTauCache> DeepTauBase::initializeGlobalCache(const edm::ParameterSet& cfg) {
    const auto graph_name_vector = cfg.getParameter<std::vector<std::string>>("graph_file");
    std::map<std::string, std::string> graph_names;
    for (const auto& entry : graph_name_vector) {
      const size_t sep_pos = entry.find(':');
      std::string entry_name, graph_file;
      if (sep_pos != std::string::npos) {
        entry_name = entry.substr(0, sep_pos);
        graph_file = entry.substr(sep_pos + 1);
      } else {
        entry_name = "";
        graph_file = entry;
      }
      graph_file = edm::FileInPath(graph_file).fullPath();
      if (graph_names.count(entry_name))
        throw cms::Exception("DeepTauCache") << "Duplicated graph entries";
      graph_names[entry_name] = graph_file;
    }
    bool mem_mapped = cfg.getParameter<bool>("mem_mapped");
    return std::make_unique<DeepTauCache>(graph_names, mem_mapped);
  }

  DeepTauCache::DeepTauCache(const std::map<std::string, std::string>& graph_names, bool mem_mapped) {
    for (const auto& graph_entry : graph_names) {
      tensorflow::SessionOptions options;
      tensorflow::setThreading(options, 1);

      const std::string& entry_name = graph_entry.first;
      const std::string& graph_file = graph_entry.second;
      if (mem_mapped) {
        memmappedEnv_[entry_name] = std::make_unique<tensorflow::MemmappedEnv>(tensorflow::Env::Default());
        const tensorflow::Status mmap_status = memmappedEnv_.at(entry_name)->InitializeFromFile(graph_file);
        if (!mmap_status.ok()) {
          throw cms::Exception("DeepTauCache: unable to initalize memmapped environment for ")
              << graph_file << ". \n"
              << mmap_status.ToString();
        }

        graphs_[entry_name] = std::make_unique<tensorflow::GraphDef>();
        const tensorflow::Status load_graph_status =
            ReadBinaryProto(memmappedEnv_.at(entry_name).get(),
                            tensorflow::MemmappedFileSystem::kMemmappedPackageDefaultGraphDef,
                            graphs_.at(entry_name).get());
        if (!load_graph_status.ok())
          throw cms::Exception("DeepTauCache: unable to load graph from ") << graph_file << ". \n"
                                                                           << load_graph_status.ToString();

        options.config.mutable_graph_options()->mutable_optimizer_options()->set_opt_level(
            ::tensorflow::OptimizerOptions::L0);
        options.env = memmappedEnv_.at(entry_name).get();

        sessions_[entry_name] = tensorflow::createSession(graphs_.at(entry_name).get(), options);

      } else {
        graphs_[entry_name].reset(tensorflow::loadGraphDef(graph_file));
        sessions_[entry_name] = tensorflow::createSession(graphs_.at(entry_name).get(), options);
      }
    }
  }

  DeepTauCache::~DeepTauCache() {
    for (auto& session_entry : sessions_)
      tensorflow::closeSession(session_entry.second);
  }

}  // namespace deep_tau

namespace deep_tau_2017 {

  float getTauID(const pat::Tau& tau, const std::string& tauID, float default_value, bool assert_input) {
    static tbb::concurrent_unordered_set<std::string> isFirstWarning;
    if (tau.isTauIDAvailable(tauID)) {
      return tau.tauID(tauID);
    } else {
      if (assert_input) {
        throw cms::Exception("DeepTauId")
            << "Exception in <getTauID>: No tauID '" << tauID << "' available in pat::Tau given as function argument.";
      }
      if (isFirstWarning.insert(tauID).second) {
        edm::LogWarning("DeepTauID") << "Warning in <getTauID>: No tauID '" << tauID
                                     << "' available in pat::Tau given as function argument."
                                     << " Using default_value = " << default_value << " instead." << std::endl;
      }
      return default_value;
    }
  }

  void MuonHitMatchV1::addMatchedMuon(const pat::Muon& muon, reco::BaseTau const& tau) {
    static constexpr int n_stations = 4;

    ++n_muons;
    const double dR2 = reco::deltaR2(tau.p4(), muon.p4());
    if (!best_matched_muon || dR2 < deltaR2_best_match) {
      best_matched_muon = &muon;
      deltaR2_best_match = dR2;
    }

    for (const auto& segment : muon.matches()) {
      if (segment.segmentMatches.empty())
        continue;
      if (n_matches.count(segment.detector()))
        ++n_matches.at(segment.detector()).at(segment.station() - 1);
    }

    if (muon.outerTrack().isNonnull()) {
      const auto& hit_pattern = muon.outerTrack()->hitPattern();
      for (int hit_index = 0; hit_index < hit_pattern.numberOfAllHits(reco::HitPattern::TRACK_HITS); ++hit_index) {
        auto hit_id = hit_pattern.getHitPattern(reco::HitPattern::TRACK_HITS, hit_index);
        if (hit_id == 0)
          break;
        if (hit_pattern.muonHitFilter(hit_id) && (hit_pattern.getHitType(hit_id) == TrackingRecHit::valid ||
                                                  hit_pattern.getHitType(hit_id == TrackingRecHit::bad))) {
          const int station = hit_pattern.getMuonStation(hit_id) - 1;
          if (station > 0 && station < n_stations) {
            std::vector<UInt_t>* muon_n_hits = nullptr;
            if (hit_pattern.muonDTHitFilter(hit_id))
              muon_n_hits = &n_hits.at(MuonSubdetId::DT);
            else if (hit_pattern.muonCSCHitFilter(hit_id))
              muon_n_hits = &n_hits.at(MuonSubdetId::CSC);
            else if (hit_pattern.muonRPCHitFilter(hit_id))
              muon_n_hits = &n_hits.at(MuonSubdetId::RPC);

            if (muon_n_hits)
              ++muon_n_hits->at(station);
          }
        }
      }
    }
  }

  unsigned MuonHitMatchV1::countMuonStationsWithMatches(size_t first_station, size_t last_station) const {
    static const std::map<int, std::vector<bool>> masks = {
        {MuonSubdetId::DT, {false, false, false, false}},
        {MuonSubdetId::CSC, {true, false, false, false}},
        {MuonSubdetId::RPC, {false, false, false, false}},
    };
    unsigned cnt = 0;
    for (unsigned n = first_station; n <= last_station; ++n) {
      for (const auto& match : n_matches) {
        if (!masks.at(match.first).at(n) && match.second.at(n) > 0)
          ++cnt;
      }
    }
    return cnt;
  }

  unsigned MuonHitMatchV1::countMuonStationsWithHits(size_t first_station, size_t last_station) const {
    static const std::map<int, std::vector<bool>> masks = {
        {MuonSubdetId::DT, {false, false, false, false}},
        {MuonSubdetId::CSC, {false, false, false, false}},
        {MuonSubdetId::RPC, {false, false, false, false}},
    };

    unsigned cnt = 0;
    for (unsigned n = first_station; n <= last_station; ++n) {
      for (const auto& hit : n_hits) {
        if (!masks.at(hit.first).at(n) && hit.second.at(n) > 0)
          ++cnt;
      }
    }
    return cnt;
  }

  template <>
  CellObjectType GetCellObjectType(const pat::Electron&) {
    return CellObjectType::Electron;
  }

  template <>
  CellObjectType GetCellObjectType(const pat::Muon&) {
    return CellObjectType::Muon;
  }

  template <>
  CellObjectType GetCellObjectType(reco::Candidate const& cand) {
    static const std::map<int, CellObjectType> obj_types = {{11, CellObjectType::PfCand_electron},
                                                            {13, CellObjectType::PfCand_muon},
                                                            {22, CellObjectType::PfCand_gamma},
                                                            {130, CellObjectType::PfCand_neutralHadron},
                                                            {211, CellObjectType::PfCand_chargedHadron}};

    auto iter = obj_types.find(std::abs(cand.pdgId()));
    if (iter == obj_types.end())
      return CellObjectType::Other;
    return iter->second;
  }

  bool CellGrid::tryGetCellIndex(double deltaEta, double deltaPhi, CellIndex& cellIndex) const {
    const auto getCellIndex = [this](double x, double maxX, double size, int& index) {
      const double absX = std::abs(x);
      if (absX > maxX)
        return false;
      double absIndex;
      if (disable_CellIndex_workaround_) {
        // CV: use consistent definition for CellIndex
        //     in DeepTauId.cc code and new DeepTau trainings
        absIndex = std::floor(absX / size + 0.5);
      } else {
        // CV: backwards compatibility with DeepTau training v2p1 used during Run 2
        absIndex = std::floor(std::abs(absX / size - 0.5));
      }
      index = static_cast<int>(std::copysign(absIndex, x));
      return true;
    };

    return getCellIndex(deltaEta, maxDeltaEta(), cellSizeEta, cellIndex.eta) &&
           getCellIndex(deltaPhi, maxDeltaPhi(), cellSizePhi, cellIndex.phi);
  }

}  // namespace deep_tau_2017

namespace deeptau_helper {
  const deep_tau::DeepTauBase::OutputCollection& GetOutputs() {
    static constexpr size_t e_index = 0, mu_index = 1, tau_index = 2, jet_index = 3;
    static const deep_tau::DeepTauBase::OutputCollection outputs_ = {
        {"VSe", deep_tau::DeepTauBase::Output({tau_index}, {e_index, tau_index})},
        {"VSmu", deep_tau::DeepTauBase::Output({tau_index}, {mu_index, tau_index})},
        {"VSjet", deep_tau::DeepTauBase::Output({tau_index}, {jet_index, tau_index})},
    };
    return outputs_;
  }

  bool isAbove(double value, double min) { return std::isnormal(value) && value > min; }

  bool calculateElectronClusterVarsV2(const pat::Electron& ele,
                                      float& cc_ele_energy,
                                      float& cc_gamma_energy,
                                      int& cc_n_gamma) {
    cc_ele_energy = cc_gamma_energy = 0;
    cc_n_gamma = 0;
    const auto& superCluster = ele.superCluster();
    if (superCluster.isNonnull() && superCluster.isAvailable() && superCluster->clusters().isNonnull() &&
        superCluster->clusters().isAvailable()) {
      for (auto iter = superCluster->clustersBegin(); iter != superCluster->clustersEnd(); ++iter) {
        const float energy = static_cast<float>((*iter)->energy());
        if (iter == superCluster->clustersBegin())
          cc_ele_energy += energy;
        else {
          cc_gamma_energy += energy;
          ++cc_n_gamma;
        }
      }
      return true;
    } else
      return false;
  }

  double getInnerSignalConeRadius(double pt) {
    static constexpr double min_pt = 30., min_radius = 0.05, cone_opening_coef = 3.;
    // This is equivalent of the original formula (std::max(std::min(0.1, 3.0/pt), 0.05)
    return std::max(cone_opening_coef / std::max(pt, min_pt), min_radius);
  }

  bool isInEcalCrack(double eta) {
    const double abs_eta = std::abs(eta);
    return abs_eta > 1.46 && abs_eta < 1.558;
  }
}  // namespace deeptau_helper
