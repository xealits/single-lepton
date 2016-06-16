runAnalysis is the chhigs/runAnalysis fully merged at 3-3-2016.
Plus the acceptance histograms.

singleLepton is the single-lepton part of it.


# Code parts


## Initialization of config parameters

**the main ones**

    // load framework libraries
    gSystem->Load ("libFWCoreFWLite");
    AutoLibraryLoader::enable ();

    const edm::ParameterSet & runProcess = edm::readPSetsFrom (argv[1])->getParameter < edm::ParameterSet > ("runProcess");

    bool debug           = runProcess.getParameter<bool>  ("debug");
    bool runSystematics  = runProcess.getParameter<bool>  ("runSystematics");
    bool saveSummaryTree = runProcess.getParameter<bool>  ("saveSummaryTree");
    bool isMC            = runProcess.getParameter<bool>  ("isMC");
    double xsec          = runProcess.getParameter<double>("xsec");
    int mctruthmode      = runProcess.getParameter<int>   ("mctruthmode");
    TString dtag         = runProcess.getParameter<std::string>("dtag");
  
    const edm::ParameterSet& myVidElectronIdConf = runProcess.getParameterSet("electronidparas");
    const edm::ParameterSet& myVidElectronMainIdWPConf = myVidElectronIdConf.getParameterSet("tight");
    const edm::ParameterSet& myVidElectronVetoIdWPConf = myVidElectronIdConf.getParameterSet("loose");

    std::vector < std::string > urls = runProcess.getUntrackedParameter < std::vector < std::string > >("input");
    TString outUrl = runProcess.getParameter<std::string>("outfile");

    lumiUtils::GoodLumiFilter goodLumiFilter(runProcess.getUntrackedParameter<std::vector<edm::LuminosityBlockRange> >("lumisToProcess", std::vector<edm::LuminosityBlockRange>()));

    if (dtag.Contains ("SingleMuon")) filterOnlySINGLEMU = true;
    ...

    bool isV0JetsMC   (isMC && (dtag.Contains ("DYJetsToLL") || dtag.Contains ("WJets")));

    // Reactivate for diboson shapes  
    // bool isMC_ZZ      (isMC && (string (dtag.Data ()).find ("TeV_ZZ") != string::npos));
    // ...

    bool isTTbarMC    (isMC && (dtag.Contains("TTJets") || dtag.Contains("_TT_") )); // Is this still useful?
    bool isPromptReco (!isMC && dtag.Contains("PromptReco")); //"False" picks up correctly the new prompt reco (2015C) and MC
    bool isRun2015B   (!isMC && dtag.Contains("Run2015B"));
    bool isNLOMC      (isMC && (dtag.Contains("amcatnlo") || dtag.Contains("powheg")) );

    //tree info
    TString dirname = runProcess.getParameter < std::string > ("dirName");

    std::vector<TString> systVars(1,"");
    if(runSystematics && isMC)
      {
      systVars.push_back("jerup" );     systVars.push_back("jerdown"   );
      systVars.push_back("jesup" );     systVars.push_back("jesdown"   );
      //systVars.push_back("lesup" );   systVars.push_back("lesdown"   );
      systVars.push_back("leffup");     systVars.push_back("leffdown"  );
      systVars.push_back("puup"  );     systVars.push_back("pudown"   );
      systVars.push_back("umetup");     systVars.push_back("umetdown" );
      systVars.push_back("btagup");     systVars.push_back("btagdown" );
      systVars.push_back("unbtagup");   systVars.push_back("unbtagdown" );
      if(isTTbarMC) {systVars.push_back("topptuncup");
         systVars.push_back("topptuncdown"); }
      //systVars.push_back(); systVars.push_back();
    
      if(isTTbarMC) { systVars.push_back("pdfup"); systVars.push_back("pdfdown"); }
      cout << "Systematics will be computed for this analysis - this will take a bit" << endl;
      }


    // TODO: what is this: allWeightsURL ... "weightsFile"??
    std::vector < std::string > allWeightsURL = runProcess.getParameter < std::vector < std::string > >("weightsFile");
    std::string weightsDir (allWeightsURL.size ()? allWeightsURL[0] : "");
    // weightsDir is not used

    //  //shape uncertainties for dibosons
    ...not used


    // ----------------------------------- jet energy scale and uncertainties 
  
    // ----------------------------------- muon energy scale and uncertainties

    // --------------------------------------- lepton efficiencies

    // --------------------------------------- b-tagging
    // --------------------------------------- electron IDs, main and veto
    // --------------------------------------- pileup weighting

    // --------------------------------------- hardcoded MET filter



## Initialization of histograms

Какие мне нужны?

Сколько всего эвентов в сете?
MC ре-формируется (ре-вейтится) чтобы соотв. данным (пайлап и др.).
Сколько эвентов в МС после ре-формировки?
(Общий фактор ещё можно корректировать.)

Далее,
channel selection & selection steps.




## Event loop, application of config parameters and selection

* weird NLO -1 weights
* pileup weight (with plus-minus)
* creeppish merging of LO and NLO sets (HT binning)
* count N good verteces
* Apply pileup reweighting
* save distributions of weights
* ?(all weight is applied -- there should be an overall integral here)
* Orthogonalize Run2015B PromptReco+17Jul15 mix
* Skip bad lumi
* apply trigger
* Apply MET filters
* load all the objects we will need to access
* "TODO: what is this??" thing
* actual particles
* merging electrons and muons
* leptons selection
  + kinematics, main and veto
  + lepton IDs and isolation
* select the taus
* JET/MET ANALYSIS
* ASSIGN CHANNEL
* Single lepton full analysis
  + Clean jet collection from selected taus
  + only selections and filling histograms




## Plotter








# Rewriting the code

## Clean-lepton, converging to Mara's config

## Event loop, application of config parameters and selection

The steps of Pietro's code with changes.

* MC shaping to data properties
  + NLO -1 weights for powheg and amcatnlo datasets -> **0 leave as is**
  + *removed* merging-stitching of LO and NLO sets (via HT/phat binning)
  + count N good verteces (used as pile-up in data?)
  + Apply pileup reweighting -> **manual reweight**
  + (TODO) HLT trigger efficiency
  + (TODO) lepton ID/Iso efficiencies
  + B-tag efficiency (needs to be checked and updated according to \url{twiki/bin/viewauth/CMS/BtagRecommendation76X})
  + (TODO?) PDF weights from \url{https://twiki.cern.ch/twiki/bin/viewauth/CMS/LHEReaderCMSSW#How_to_use_weights}
  + save distributions of weights -> **1**
* basic event selection
  + remove Run2015B[Orthogonalize Run2015B PromptReco+17Jul15 mix] **0**
  + Skip bad lumi -> **check for lumicert for new datasets 1**
    using: `Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_v2.txt`
    (RunD is about 2160.8 pb^-1 RunC ~ 17.2 pb^-1)
    (old: `Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt`)
  + apply trigger -> **new triggers**
    muons --- `HLT_IsoMu20` or `HLT_IsoTkMu20` for data and MC
    (to update to `HLT_IsoMu18)
    electrons --- `HLT_Ele23_WPLoose_Gsf` for data and MC
    The dilepton HLTs are not used --- should consider adding them (people regard this option as efficiency improvement \url{https://hypernews.cern.ch/HyperNews/CMS/get/top-trigger/180.html})
  + *Apply MET filters*
    using now DoubleEG_csc2015.txt, MuonEG_csc2015.txt etc
    TODO: check "MET-filter triggers" for 76X in https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2015#ETmiss_filters
* load all the objects we will need to access
* "TODO: what is this??" thing -> **0 commented out**
   (it should be the electron-muon split)
* actual particles:
  - muons,
  - electrons,
  - jets,
  - METs (collection of MET -- there are different MET algorithms!)
     **we use the 0th MET of slimmedMETs collection in 76x MINIAODs v2 -- it is type 1 MET**
  - taus
* not merging electrons and muons
* leptons selection
  + muon corrections are applied with rochester correction procedure
    implemented in https://github.com/cms2l2v/2l2v_fwk/blob/master/interface/rochcor2015.h
  + electron corrections applied with CMSSW
    `EgammaAnalysis/ElectronTools/interface/ElectronEnergyCalibratorRun2.h`
    https://github.com/cms-sw/cmssw/blob/CMSSW_8_1_X/EgammaAnalysis/ElectronTools/interface/ElectronEnergyCalibratorRun2.h
    according to https://twiki.cern.ch/twiki/bin/viewauth/CMS/EGMSmearer
  + kinematics, good and veto, lepton IDs and isolation -> **new isolation**
    - muons:
      good: P_T > 26, eta < 2.4, tight muon
      veto: P_T > 10, eta < 2.5, loose muon

      IDs are according to
      https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
      These IDs are assigned to muons in data for CMSSW_7_4_2 and above.
      They are combinations of selectiors:
        - loose muon:
          recoMu.isPFMuon() & (recoMu.isGlobalMuon() || recoMu.isTrackerMuon())
        - tight muon:
          recoMu.isGlobalMuon() &
          recoMu.isPFMuon() &
          recoMu.globalTrack()->normalizedChi2() < 10. &
          recoMu.globalTrack()->hitPattern().numberOfValidMuonHits() > 0 &
          recoMu.numberOfMatchedStations() > 1 &
          fabs(recoMu.muonBestTrack()->dxy(vertex->position())) < 0.2
          Or dB() < 0.2 on pat::Muon [1]
          fabs(recoMu.muonBestTrack()->dz(vertex->position())) < 0.5
          recoMu.innerTrack()->hitPattern().numberOfValidPixelHits() > 0
          recoMu.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5

      Muons isolation is done according to
      https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Muon_Isolation
      by the scheme called
      "PF-based combined relative isolation with deltaBeta correction."

    - electrons:
      good: P_T > 30, eta < 2.4,
      veto: P_T > 15, eta < 2.5,
      both have window at leta > 1.4442 && leta < 1.5660
      IDs and isolation are done with cut-based procedure according to
      https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Spring15_selection_25ns
      (Isolation twiki requested)

* select the taus -> **0 leave as is**
  + 4 tau ID discriminators (**check if up to date**) from https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV:
    decayModeFinding > 0.5
    byMediumCombinedIsolationDeltaBetaCorr3Hits > 0.5
    againstMuonTight3 > 0.5
    (old againstElectronMediumMVA6 > 0.5)
    againstElectronTightMVA6 > 0.5
  + pixel hits cut (*"should be available out of the mox in new MINIAOD"* --- ?):
    basically now we check if tau.signalChargedHadrCands
    have at least 1 element with numberOfPixelHits > 0
  + tau pt > 20, eta < 2.3
  + overlap with selLeptons $\delta R > 0.4 $
* jets
  + corrections `Fall15_25nsV2_*`
    - uncorrect with `jet.correctedP4("Uncorrected")`
    - apply `FactorizedJetCorrector` for JES (jet energy scale) corrections
      https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetEnCorFWLite
    - smear JER (jet energy resolution) in MC
      https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#JER_Scaling_factors_and_Uncertai
  + individual jet selection
    - IDs from https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_data
    loose jet is required
    - pt > 30, eta < 2.5
  + dphijmet = fabs( deltaPhi(curr_jet, met) ) -- and save the min
  + cross-clean of leptons and taus with deltaR > 0.4
* b-jets via "pfCombinedInclusiveSecondaryVertexV2BJetTags" > 0.89 (medium working point)
  from https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation74X
  (move to https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation76X)

* MET
  + MET 0 from slimmedMETs is used (TODO: confirm correction type)
    https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2015#ETmiss
  + (disabled new procedure due to large mismatch) propagate lepton corrections to MET
  + propagate jet and lepton corrections to MET with getMETvariations
    (TODO: check getMETvariations -- there is some electron-muon difference)
    (it seems getMETvariations does not propagate the jet corrections,
     it only variates the MET for systematics)
* ASSIGN CHANNEL
  + single lepton:
    only 1 good lepton, no veto, number of good PV > 0
    2 or more jets
    MET $p_T > 40 GeV$
    1 or more b-taged jets
    Only 1 tau
    Oposite electric charge of tau and lepton

  + double lepton:
    Only 2 good leptons, no veto, number of good PV > 0
    Mass of dilepton system $> 12 GeV$ and Z-window of $15 GeV$ for $ee$, $\mu\mu$ around $91 GeV$
    2 or more jets
    MET $p_T > 40 GeV$
    Oposite sign of leptons
    1 or more b-jets

* Single lepton full analysis
  + *Clean jet collection from selected taus* moved it up to common selection
  + only selections and filling histograms
  + our selection -> **0 + Mara's selection: 1 lepton, 4 jets, 2btags**

-- more steps?

Other Mara's steps:

* MC normalization
  + MC weights twiki/bin/viewauth/CMS/LHEReaderCMSSW#How_to_use_weights (PDF and others?)
  + pile-up weighting -- **the same now**
  + Muon eff, isolation, ID, trigger (TODO)
  + b-tagging efficiencies twiki/bin/viewauth/CMS/BtagRecommendation76X (TO UPDATE)
* different datasets (**using them now**)

What about trigger efficiency? Do data and MC really match above the threshold?

Are electron and muon datasets of the same run orthogonal?

