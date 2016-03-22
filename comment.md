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






