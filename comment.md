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
      ...
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

* weird NLO -1 weights -> **0 leave as is**
* pileup weight (with plus-minus) -> *manual weights Y*
* *removed*[creeppish merging of LO and NLO sets (HT binning)]
* count N good verteces
* Apply pileup reweighting -> *manual reweight Y*
* save distributions of weights -> **1**
* ?(all weight is applied -- there should be an overall integral here)
* remove Run2015B[Orthogonalize Run2015B PromptReco+17Jul15 mix] **0**
* Skip bad lumi -> **check for lumicert for new datasets 1**
* apply trigger -> new triggers **1**
* Apply MET filters -> **bug in metFilter, running bug-less * cleanLepton**
* load all the objects we will need to access
* "TODO: what is this??" thing -> **0 commented out**
* actual particles:
  - muons,
  - electrons,
  - jets,
  - gammas,
  - METs (collection of MET -- there are different MET algorithms!)
     **check which MET we use**
  - taus
* merging electrons and muons
* leptons selection
  + apply muon corrections, muCor
  + kinematics, main and veto -> new threshold **1**
  + lepton IDs and isolation -> new isolation **1**
* select the taus -> **0 ?? leave as is**
* JET/MET ANALYSIS -> **0 ?? leave as is** -> update b-jets
  + it only selects jets, some of them -- as b-tagged
  + only 1 parameter is obtained from MET
  + updateJEC
  + newMet = getMETvariations
  + selecting jets:
      pt, eta, **mctruth (?)**, cross-clean with l/gamma dR, jet ID
  + dphijmet = fabs( deltaPhi(curr_jet, met) ) -- and save the min
  + b-tagging via hasCSVtag
* ASSIGN CHANNEL
* Single lepton full analysis
  + Clean jet collection from selected taus
  + only selections and filling histograms
  + 6 selection steps -> **0 + Mara's selection: 1 lepton, 4 jets, 2btags**

-- there were more steps, smearing muon momentum,
steps in taus and jets.

Other Mara's steps:

* in MC normalization
  + MC weights twiki/bin/viewauth/CMS/LHEReaderCMSSW#How_to_use_weights
  + pile-up weighting -- **the same now**
  + Muon eff, isolation, ID, trigger (??)
  + b-tagging efficiencies twiki/bin/viewauth/CMS/BtagRecommendation76X
* different datasets


