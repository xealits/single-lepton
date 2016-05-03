//
// Pietro Vischia, <pietro.vischia@gmail.com>
//
// ttbar and charged Higgs analyses
//

#include <iostream>
#include <boost/shared_ptr.hpp>

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

//Load here all the dataformat that we will need
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "GeneratorInterface/LHEInterface/interface/LHEEvent.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
//#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h" //for svfit

#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/HiggsUtils.h"
#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/llvv_fwk/interface/TMVAUtils.h"
#include "UserCode/llvv_fwk/interface/LeptonEfficiencySF.h"
#include "UserCode/llvv_fwk/interface/PDFInfo.h"
#include "UserCode/llvv_fwk/interface/MuScleFitCorrector.h"
#include "UserCode/llvv_fwk/interface/BtagUncertaintyComputer.h"
#include "UserCode/llvv_fwk/interface/GammaWeightsHandler.h"

#include "UserCode/llvv_fwk/interface/PatUtils.h"


#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TEventList.h"
#include "TROOT.h"
#include "TNtuple.h"
#include <Math/VectorUtil.h>

using namespace std;

namespace utils
	{
	namespace cmssw
		{
		std::vector<double> smearJER(double pt, double eta, double genPt)
			{
			std::vector<double> toReturn(3,pt);
			if(genPt<=0) return toReturn;
			
			// FIXME: These are the 8 TeV values.
			//
			eta=fabs(eta);
			double ptSF(1.0), ptSF_err(0.06);
			if(eta<0.5)                  { ptSF=1.052; ptSF_err=sqrt(pow(0.012,2)+pow(0.5*(0.062+0.061),2)); }
			else if(eta>=0.5 && eta<1.1) { ptSF=1.057; ptSF_err=sqrt(pow(0.012,2)+pow(0.5*(0.056+0.055),2)); }
			else if(eta>=1.1 && eta<1.7) { ptSF=1.096; ptSF_err=sqrt(pow(0.017,2)+pow(0.5*(0.063+0.062),2)); }
			else if(eta>=1.7 && eta<2.3) { ptSF=1.134; ptSF_err=sqrt(pow(0.035,2)+pow(0.5*(0.087+0.085),2)); }
			else if(eta>=2.3 && eta<5.0) { ptSF=1.288; ptSF_err=sqrt(pow(0.127,2)+pow(0.5*(0.155+0.153),2)); }
			
			toReturn[0]=TMath::Max(0.,(genPt+ptSF*(pt-genPt)));
			toReturn[1]=TMath::Max(0.,(genPt+(ptSF+ptSF_err)*(pt-genPt)));
			toReturn[2]=TMath::Max(0.,(genPt+(ptSF-ptSF_err)*(pt-genPt)));
			return toReturn;
			}

		//
		std::vector<double> smearJES(double pt, double eta, JetCorrectionUncertainty *jecUnc)
			{
			jecUnc->setJetEta(eta);
			jecUnc->setJetPt(pt);
			double relShift=fabs(jecUnc->getUncertainty(true));
			std::vector<double> toRet;
			toRet.push_back((1.0+relShift)*pt);
			toRet.push_back((1.0-relShift)*pt);
			return toRet;
			}
		
		void updateJEC(pat::JetCollection &jets, FactorizedJetCorrector *jesCor, JetCorrectionUncertainty *totalJESUnc, float rho, int nvtx,bool isMC)
			{
			for(size_t ijet=0; ijet<jets.size(); ijet++)
				{
				pat::Jet jet = jets[ijet];
				//correct JES
				LorentzVector rawJet = jet.correctedP4("Uncorrected");
				//double toRawSF=jet.correctedJet("Uncorrected").pt()/jet.pt();
				//LorentzVector rawJet(jet*toRawSF);
				jesCor->setJetEta(rawJet.eta());
				jesCor->setJetPt(rawJet.pt());
				jesCor->setJetA(jet.jetArea());
				jesCor->setRho(rho);
				jesCor->setNPV(nvtx);
				double newJECSF=jesCor->getCorrection();
				rawJet *= newJECSF;
				//jet.SetPxPyPzE(rawJet.px(),rawJet.py(),rawJet.pz(),rawJet.energy());
				jet.setP4(rawJet);

				//smear JER
				double newJERSF(1.0);
				if(isMC)
					{
					const reco::GenJet* genJet=jet.genJet();
					double genjetpt( genJet ? genJet->pt(): 0.);
					std::vector<double> smearJER=utils::cmssw::smearJER(jet.pt(),jet.eta(),genjetpt);
					newJERSF=smearJER[0]/jet.pt();
					rawJet *= newJERSF;
					//jet.SetPxPyPzE(rawJet.px(),rawJet.py(),rawJet.pz(),rawJet.energy());
					jet.setP4(rawJet);
					// FIXME: change the way this is stored (to not storing it)
					// //set the JER up/down alternatives 
					// jets[ijet].setVal("jerup",   smearJER[1] );
					// jets[ijet].setVal("jerdown", smearJER[2] );
					}

				// FIXME: change the way this is stored (to not storing it)
				////set the JES up/down pT alternatives
				//std::vector<float> ptUnc=utils::cmssw::smearJES(jet.pt(),jet.eta(), totalJESUnc);
				//jets[ijet].setVal("jesup",    ptUnc[0] );
				//jets[ijet].setVal("jesdown",  ptUnc[1] );

				// FIXME: this is not to be re-set. Check that this is a desired non-feature.
				// i.e. check that the uncorrectedJet remains the same even when the corrected momentum is changed by this routine. 
				//to get the raw jet again
				//jets[ijet].setVal("torawsf",1./(newJECSF*newJERSF));  
				}
			}

		enum METvariations { NOMINAL, JERUP, JERDOWN, JESUP, JESDOWN, UMETUP, UMETDOWN, LESUP, LESDOWN };    

		std::vector<LorentzVector> getMETvariations(LorentzVector &rawMETP4, pat::JetCollection &jets, std::vector<patUtils::GenericLepton> &leptons,bool isMC)
			{
			std::vector<LorentzVector> newMetsP4(9,rawMETP4);
			if(!isMC) return newMetsP4;
			
			LorentzVector nullP4(0,0,0,0);
			//recompute the clustered and unclustered fluxes with energy variations
			for(size_t ivar=1; ivar<=8; ivar++)
				{
				//leptonic flux
				LorentzVector leptonFlux(nullP4), lepDiff(nullP4);
				for(size_t ilep=0; ilep<leptons.size(); ilep++)
					{
					LorentzVector lepton = leptons[ilep].p4();
					double varSign( (ivar==LESUP ? 1.0 : (ivar==LESDOWN ? -1.0 : 0.0) ) );
					int id( abs(leptons[ilep].pdgId()) );
					double sf(1.0);
					if(id==13) sf=(1.0+varSign*0.01);
					if(id==11)
						{
						if(fabs(leptons[ilep].eta())<1.442) sf=(1.0+varSign*0.02);
						else                                sf=(1.0-varSign*0.05);
						}
					leptonFlux += lepton;
					lepDiff += (sf-1)*lepton;
					}
			
				//clustered flux
				LorentzVector jetDiff(nullP4), clusteredFlux(nullP4);
				for(size_t ijet=0; ijet<jets.size(); ijet++)
					{
					if(jets[ijet].pt()==0) continue;
					double jetsf(1.0);
					// FIXME: change the way this is stored (to not storing it)              
					/// if(ivar==JERUP)   jetsf=jets[ijet].getVal("jerup")/jets[ijet].pt();
					/// if(ivar==JERDOWN) jetsf=jets[ijet].getVal("jerdown")/jets[ijet].pt();
					/// if(ivar==JESUP)   jetsf=jets[ijet].getVal("jesup")/jets[ijet].pt();
					/// if(ivar==JESDOWN) jetsf=jets[ijet].getVal("jesdown")/jets[ijet].pt();
					//LorentzVector newJet( jets[ijet] ); newJet *= jetsf;
					LorentzVector newJet = jets[ijet].p4(); newJet *= jetsf;
					jetDiff       += (newJet-jets[ijet].p4());
					clusteredFlux += jets[ijet].p4();
					}
				LorentzVector iMet=rawMETP4-jetDiff-lepDiff;

				//unclustered flux
				if(ivar==UMETUP || ivar==UMETDOWN)
					{
					LorentzVector unclusteredFlux=-(iMet+clusteredFlux+leptonFlux);
					unclusteredFlux *= (ivar==UMETUP ? 1.1 : 0.9); 
					iMet = -clusteredFlux -leptonFlux - unclusteredFlux;
					}
			
				//save new met
				newMetsP4[ivar]=iMet;
				}

			//all done here
			return newMetsP4;
			}
		}
	}

bool passPFJetID(std::string label, pat::Jet jet)
{
	bool passID(false); 
	
	float rawJetEn(jet.correctedJet("Uncorrected").energy() );

	double eta=jet.eta();
 
	float nhf( (jet.neutralHadronEnergy() + jet.HFHadronEnergy())/rawJetEn );
	float nef( jet.neutralEmEnergy()/rawJetEn );
	float cef( jet.chargedEmEnergy()/rawJetEn );
	float chf( jet.chargedHadronEnergy()/rawJetEn );
	float nch    = jet.chargedMultiplicity();
	float nconst = jet.chargedMultiplicity()+jet.neutralMultiplicity();
	float muf(jet.muonEnergy()/rawJetEn); 

	// Set of cuts from the POG group: https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_data
	if (label=="Loose")
		passID = ( ((nhf<0.99 && nef<0.99 && nconst>1) && ((abs(eta)<=2.4 && chf>0 && nch>0 && cef<0.99) || abs(eta)>2.4)) && abs(eta)<=3.0 );
	if (label=="Tight")
		passID = ( ((nhf<0.90 && nef<0.90 && nconst>1) && ((abs(eta)<=2.4 && chf>0 && nch>0 && cef<0.90) || abs(eta)>2.4)) && abs(eta) <=3.0);
	
	// Should be added the abs(eta)>3.0 part, but we never consider such jets, so... Meh!
	
	return passID; 
}



bool hasLeptonAsDaughter(const reco::GenParticle p)
{
	bool foundL(false);
	if(p.numberOfDaughters()==0) return foundL;

	// cout << "Particle " << p.pdgId() << " with status " << p.status() << " and " << p.numberOfDaughters() << endl;
	const reco::Candidate *part = &p;
	// loop on the daughter particles to check if it has an e/mu as daughter
	while ((part->numberOfDaughters()>0))
		{
		const reco::Candidate* DaughterPart = part->daughter(0);
		// cout << "\t\t Daughter: " << DaughterPart->pdgId() << " with status " << DaughterPart->status() << endl;
		if (fabs(DaughterPart->pdgId()) == 11 || fabs(DaughterPart->pdgId() == 13))
			{
			foundL = true;
			break;
			}
		part=DaughterPart;
		}
	return foundL;
}


bool hasWasMother(const reco::GenParticle  p)
{
	bool foundW(false);
	if(p.numberOfMothers()==0) return foundW;
	const reco::Candidate* part =&p; // (p.mother());
	// loop on the mother particles to check if it has a W as mother
	while ((part->numberOfMothers()>0))
		{
		const reco::Candidate* MomPart =part->mother();
		if (fabs(MomPart->pdgId())==24)
			{
			foundW = true;
			break;
			}
		part = MomPart;
		}
	return foundW;
}

bool hasTauAsMother(const reco::GenParticle  p)
{
	bool foundTau(false);
	if (p.numberOfMothers()==0) return foundTau;
	const reco::Candidate* part = &p; //(p.mother());
	// loop on the mother particles to check if it has a tau as mother
	while ((part->numberOfMothers()>0))
		{
		const reco::Candidate* MomPart =part->mother();
		if (fabs(MomPart->pdgId())==15)// && MomPart->status() == 2) // Not sure the status check is needed.
			{
			foundTau = true;
			break;
			}
		part = MomPart;
		}
	return foundTau;
}

















int main (int argc, char *argv[])
{
//##############################################
//########2    GLOBAL INITIALIZATION     ########
//##############################################

// check arguments
if (argc < 2)
	{
		std::cout << "Usage : " << argv[0] << " parameters_cfg.py" << std::endl;
		exit (0);
	}
	
// load framework libraries
gSystem->Load ("libFWCoreFWLite");
AutoLibraryLoader::enable ();
	
// configure the process
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
	
VersionedPatElectronSelector electronVidMainId(myVidElectronMainIdWPConf);
VersionedPatElectronSelector electronVidVetoId(myVidElectronVetoIdWPConf);
	
TString suffix = runProcess.getParameter < std::string > ("suffix");
std::vector < std::string > urls = runProcess.getUntrackedParameter < std::vector < std::string > >("input");
//TString baseDir = runProcess.getParameter < std::string > ("dirName");
//  if (mctruthmode != 0) //FIXME
//    {
//      outFileUrl += "_filt";
//      outFileUrl += mctruthmode;
//    }
TString outUrl = runProcess.getParameter<std::string>("outfile");
	
// Good lumi mask
lumiUtils::GoodLumiFilter goodLumiFilter(runProcess.getUntrackedParameter<std::vector<edm::LuminosityBlockRange> >("lumisToProcess", std::vector<edm::LuminosityBlockRange>()));

bool
	filterOnlySINGLEE  (false),
	filterOnlySINGLEMU (false);
if (!isMC)
	{
	if (dtag.Contains ("SingleMuon")) filterOnlySINGLEMU = true;
	if (dtag.Contains ("SingleEle"))  filterOnlySINGLEE  = true;
	}

bool isV0JetsMC   (isMC && (dtag.Contains ("DYJetsToLL") || dtag.Contains ("WJets")));
// Reactivate for diboson shapes  
// bool isMC_ZZ      (isMC && (string (dtag.Data ()).find ("TeV_ZZ") != string::npos));
// bool isMC_WZ      (isMC && (string (dtag.Data ()).find ("TeV_WZ") != string::npos));
	
bool isTTbarMC    (isMC && (dtag.Contains("TTJets") || dtag.Contains("_TT_") )); // Is this still useful?
bool isPromptReco (!isMC && dtag.Contains("PromptReco")); //"False" picks up correctly the new prompt reco (2015C) and MC
bool isRun2015B   (!isMC && dtag.Contains("Run2015B"));
bool isNLOMC      (isMC && (dtag.Contains("amcatnlo") || dtag.Contains("powheg")) );
	
TString outTxtUrl = outUrl + ".txt";
FILE *outTxtFile = NULL;
if (!isMC) outTxtFile = fopen (outTxtUrl.Data (), "w");
printf ("TextFile URL = %s\n", outTxtUrl.Data ());

//tree info
TString dirname = runProcess.getParameter < std::string > ("dirName");

//systematics
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
	if(isTTbarMC) {systVars.push_back("topptuncup"); systVars.push_back("topptuncdown"); }
	//systVars.push_back(); systVars.push_back();

	if(isTTbarMC) { systVars.push_back("pdfup"); systVars.push_back("pdfdown"); }
	cout << "Systematics will be computed for this analysis - this will take a bit" << endl;
	}

size_t nSystVars(systVars.size());
	

// TODO: what is this: allWeightsURL ... "weightsFile"??
std::vector < std::string > allWeightsURL = runProcess.getParameter < std::vector < std::string > >("weightsFile");
std::string weightsDir (allWeightsURL.size ()? allWeightsURL[0] : "");
// weightsDir is not used
//  //shape uncertainties for dibosons
//  std::vector<TGraph *> vvShapeUnc;
//  if(isMC_ZZ || isMC_WZ)
//    {
//      TString weightsFile=weightsDir+"/zzQ2unc.root";
//      TString dist("zzpt");
//      if(isMC_WZ) { weightsFile.ReplaceAll("zzQ2","wzQ2"); dist.ReplaceAll("zzpt","wzpt"); }
//      gSystem->ExpandPathName(weightsFile);
//      TFile *q2UncF=TFile::Open(weightsFile);
//      if(q2UncF){
//    vvShapeUnc.push_back( new TGraph( (TH1 *)q2UncF->Get(dist+"_up") ) );
//    vvShapeUnc.push_back( new TGraph( (TH1 *)q2UncF->Get(dist+"_down") ) );
//    q2UncF->Close();
//      }
//    }




	
//##############################################
//######## GET READY FOR THE EVENT LOOP ########
//##############################################
size_t totalEntries(0);

TFile* summaryFile = NULL;
TTree* summaryTree = NULL; //ev->;

//  
//  if(saveSummaryTree)
//    {
//      TDirectory* cwd = gDirectory;
//      std::string summaryFileName(outUrl); 
//      summaryFileName.replace(summaryFileName.find(".root", 0), 5, "_summary.root");
//      
//      summaryFile = new TFile(summaryFileName.c_str() "recreate");
//      
//      summaryTree = new TTree("Events", "Events");
//      KEY: TTreeMetaData;1
//      KEY: TTreeParameterSets;1
//      KEY: TTreeParentage;1
//      KEY: TTreeEvents;1
//      KEY: TTreeLuminosityBlocks;1
//      KEY: TTreeRuns;
//      summaryTree->SetDirectory(summaryFile);  // This line is probably not needed
//      
//      summmaryTree->Branch(
//
//      cwd->cd();
//    }
//


//MC normalization (to 1/pb)
if(debug) cout << "DEBUG: xsec: " << xsec << endl;

// ------------------------------------- jet energy scale and uncertainties 
TString jecDir = runProcess.getParameter < std::string > ("jecDir");
gSystem->ExpandPathName (jecDir);
FactorizedJetCorrector *jesCor = utils::cmssw::getJetCorrector (jecDir, isMC);
JetCorrectionUncertainty *totalJESUnc = new JetCorrectionUncertainty ((jecDir + "/MC_Uncertainty_AK4PFchs.txt").Data ());
	
// ------------------------------------- muon energy scale and uncertainties
MuScleFitCorrector *muCor = NULL;
// FIXME: MuScle fit corrections for 13 TeV not available yet (more Zs are needed) getMuonCorrector (jecDir, dtag);

// --------------------------------------- lepton efficiencies
LeptonEfficiencySF lepEff;
	
// --------------------------------------- b-tagging 
// TODO: move all these numbers to where they are applied??
// btagMedium is used twice in the code
// merge those tagging procedures and eliminated the variable?

// Prescriptions taken from: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation74X

// b-tagging working points for 50ns 
//   (pfC|c)ombinedInclusiveSecondaryVertexV2BJetTags
//      v2CSVv2L 0.605
//      v2CSVv2M 0.890
//      v2CSVv2T 0.970
double
	btagLoose(0.605), // not used anywhere in the code
	btagMedium(0.890), // used twice in the code
	btagTight(0.970); // not used anywhere in the code

//b-tagging: scale factors
//beff and leff must be derived from the MC sample using the discriminator vs flavor
//the scale factors are taken as average numbers from the pT dependent curves see:
//https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG#2012_Data_and_MC_EPS13_prescript
BTagSFUtil btsfutil;
double beff(0.68), sfb(0.99), sfbunc(0.015);
double leff(0.13), sfl(1.05), sflunc(0.12);

// Btag SF and eff from https://indico.cern.ch/event/437675/#preview:1629681
sfb = 0.861;
// sbbunc =;
beff = 0.559;


// ------------------------------ electron IDs
// does not appear anywhere in the code at all
// TString
	// electronIdMainTag("cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
	// electronIdVetoTag("cutBasedElectronID-Spring15-25ns-V1-standalone-tight");

// --------------------------------------- pileup weighting
// pile-up is done directly with direct_pileup_reweight
std::vector<double> direct_pileup_reweight = runProcess.getParameter < std::vector < double >>("pileup_reweight_direct");
	
gROOT->cd ();                 //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE

// --------------------------------------- hardcoded MET filter
patUtils::MetFilter metFiler;
if(!isMC)
	{
	metFiler.FillBadEvents(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/MetFilter/DoubleEG_RunD/DoubleEG_csc2015.txt");
	metFiler.FillBadEvents(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/MetFilter/DoubleEG_RunD/DoubleEG_ecalscn1043093.txt");
	metFiler.FillBadEvents(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/MetFilter/DoubleMuon_RunD/DoubleMuon_csc2015.txt");
	metFiler.FillBadEvents(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/MetFilter/DoubleMuon_RunD/DoubleMuon_ecalscn1043093.txt");
	metFiler.FillBadEvents(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/MetFilter/MuonEG_RunD/MuonEG_csc2015.txt");
	metFiler.FillBadEvents(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/MetFilter/MuonEG_RunD/MuonEG_ecalscn1043093.txt");
	}


// ----------------------------
// So here we got all the parameters from the config




//##############################################
//########    INITIATING HISTOGRAMS     ########
//##############################################

// Removed the SmartSelectionMonitor
// SmartSelectionMonitor mon;

TH1D* singlelep_ttbar_initialevents  = (TH1D*) new TH1D("singlelep_ttbar_init",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  ); 
TH1D* singlelep_ttbar_preselectedevents = (TH1D*) new TH1D("singlelep_ttbar_presele",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  ); 
TH1D* singlelep_ttbar_selected_mu_events = (TH1D*) new TH1D("singlelep_ttbar_sele_mu",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  ); 
TH1D* singlelep_ttbar_selected_el_events = (TH1D*) new TH1D("singlelep_ttbar_sele_el",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  ); 
TH1D* singlelep_ttbar_selected2_mu_events = (TH1D*) new TH1D("singlelep_ttbar_sele2_mu",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  ); 
TH1D* singlelep_ttbar_selected2_el_events = (TH1D*) new TH1D("singlelep_ttbar_sele2_el",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  ); 

TH1D* singlelep_ttbar_maraselected_mu_events = (TH1D*) new TH1D("singlelep_ttbar_marasele_mu",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  ); 
TH1D* singlelep_ttbar_maraselected_el_events = (TH1D*) new TH1D("singlelep_ttbar_marasele_el",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  ); 


// Kinematic parameters of the decay
TLorentzVector pl, plb, pb, pbb, prest;

// -------------------------------
// Here the output histograms and other object should be initialized




//##############################################
//########           EVENT LOOP         ########
//##############################################
//loop on all the events
printf ("Progressing Bar     :0%%       20%%       40%%       60%%       80%%       100%%\n");

int nMultiChannel(0);
FILE *csv_out;
string FileName = ((outUrl.ReplaceAll(".root",""))+".csv").Data();
csv_out = fopen(FileName.c_str(), "w");

fprintf(csv_out, "Headers\n");


for(size_t f=0; f<urls.size();++f){
	fprintf(csv_out, "Processing file: %s\n", urls[f].c_str());
	TFile* file = TFile::Open(urls[f].c_str());
	fwlite::Event ev(file);
	printf ("Scanning the ntuple %2lu/%2lu :",f+1, urls.size());


	unsigned int n_good_muon_pairs = 0;
	unsigned int n_good_electron_pairs = 0;
	unsigned int muon_selection_control[11]     = {0,0,0,0,0, 0,0,0,0,0, 0};
	unsigned int electron_selection_control[11] = {0,0,0,0,0, 0,0,0,0,0, 0};

	// acceptance parameters
	int iev(0); // number of events
	double sum_weights_raw = 0; // sum of raw weights
	double sum_weights = 0; // sum of final weights
	double crossel_sum_weights_raw = 0; // crossel
	double crossel_sum_weights = 0;
	double oursel_sum_weights_raw = 0; // oursel
	double oursel_sum_weights = 0;
	double marasel_sum_weights_raw = 0; // marasel
	double marasel_sum_weights = 0;
	unsigned int negative_event_nvtx[100];
	unsigned int positive_event_nvtx[100];
	double negative_event_pernvtx_weight[100];
	double positive_event_pernvtx_weight[100];
	for (int i=0; i<100; i++)
		{
		negative_event_nvtx[i] = 0;
		positive_event_nvtx[i] = 0;
		negative_event_pernvtx_weight[i] = 0;
		positive_event_pernvtx_weight[i] = 0;
		}

	int treeStep (ev.size()/50);

	for (ev.toBegin(); !ev.atEnd(); ++ev)
		{
		singlelep_ttbar_initialevents->Fill(1);
		iev++;
		totalEntries++;
		if (iev % treeStep == 0)
			{
			printf (".");
			if(!debug) fflush (stdout); // Otherwise debug messages are flushed
			}

		edm::EventBase const & myEvent = ev;

		// ---------------------------------- MC shaping to data
		// MC is weighted according to distributions of a bunch of data properties

		// NLO -1 corrections
		double weightGen(1.);

		// there are also the (dissabled now, since NLO samples are used) HT-binned and pthat-binned stitching of LO and NLO

		// ---------------------------------- pileup weight
		double puWeight         (1.0);

		// rawWeight is everything but Pile-Up
		double rawWeight        (1.0);

		// final weight of the event
		double weight           (1.0);
		// and systematic corrections? TODO: check how TotalWeight_plus is used?
		double TotalWeight_plus (1.0);
		double TotalWeight_minus(1.0);


		// ---------------------------------- these are weird NLO -1 weights
		// TODO: figure out how exactly they correct for NLO
		// Take into account the negative weights from some NLO generators (otherwise some phase space will be double counted)
		if(isNLOMC)
			{

			fwlite::Handle<GenEventInfoProduct> evt;
			evt.getByLabel(ev, "generator");
			if(evt.isValid())
				{
				weightGen = (evt->weight() > 0 ) ? 1. : -1. ;
				}

			// FIXME: this is for PDF uncertainties, must reactivate it at some point.
			//fwlite::Handle<LHEEventProduct> lheEvtProd;
			//lheEvtProd.getByLabel(ev, "externalLHEProducer");
			//if(lheEvtProd.isValid())
			//  {
			//    weightLhe=lheEvtProd->originalXWGTUP();
			//    
			//   //for(unsigned int i=0; i<evet->weights().size();i++){
			//   //  double asdde=evet->weights()[i].wgt;
			//   //  EventInfo.ttbar_w[EventInfo.ttbar_nw]=EventInfo.ttbar_w[0]*asdde/asdd;
			//   //  EventInfo.ttbar_nw++;
			//   //}
			//  }
			//cout << "Event " << iev << " has genweight: " << weightGen << " and LHE weight " << weightLhe << endl;

			}


		std::vector < TString > tags (1, "all"); // Inclusive inclusiveness

		//
		// DERIVE WEIGHTS TO APPLY TO SAMPLE
		//


		// -------------- it should be the creeppish merging of LO and NLO sets
		// not used now at all?
		// This must remain deactivated if you use HT-binned samples (it was for pthat-binned samples)
		// if (isV0JetsMC)
		//   {
		//   fwlite::Handle < LHEEventProduct > lheEPHandle;
		//   lheEPHandle.getByLabel (ev, "externalLHEProducer");
		//   mon.fillHisto ("nup", "", lheEPHandle->hepeup ().NUP, 1);
		//   if (lheEPHandle->hepeup ().NUP > 5)  continue;
		//   mon.fillHisto ("nupfilt", "", lheEPHandle->hepeup ().NUP, 1);
		//   }

		// HT-binned samples stitching: https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorking2015#MC_and_data_samples
		/*
		if(isV0JetsMC)
			{
			// access generator level HT               
			fwlite::Handle<LHEEventProduct> lheEventProduct;
			lheEventProduct.getByLabel(ev, "externalLHEProducer");
			//edm::Handle<LHEEventProduct> lheEventProduct;
			//ev.getByLabel( 'externalLHEProducer', lheEventProduct);
			const lhef::HEPEUP& lheEvent = lheEventProduct->hepeup(); 
			std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;
			double lheHt = 0.;
			size_t numParticles = lheParticles.size();
			for ( size_t idxParticle = 0; idxParticle < numParticles; ++idxParticle )
				{
				int absPdgId = TMath::Abs(lheEvent.IDUP[idxParticle]);
				int status = lheEvent.ISTUP[idxParticle];
				if ( status == 1 && ((absPdgId >= 1 && absPdgId <= 6) || absPdgId == 21) )
					{
					// quarks and gluons
					lheHt += TMath::Sqrt(TMath::Power(lheParticles[idxParticle][0], 2.) + TMath::Power(lheParticles[idxParticle][1], 2.));
					// first entry is px, second py
					}                                        
				}
			if(debug) cout << "Sample: " << dtag << ", lheHt: " << lheHt << ", scale factor from spreadsheet: " << patUtils::getHTScaleFactor(dtag, lheHt) << endl;
			// getHTScaleFactor works on combining several LO datasets with NLO
			// now one 1 NLO dataset is used for both WJets and DYJets
			// thus it is commented out here
			//weightGen *=   patUtils::getHTScaleFactor(dtag, lheHt);
			}
		*/

		weight *= weightGen;
		rawWeight *=weightGen;
				
		// ------------------------------- count N good verteces
		// needed for particle selection/event classification later
		// and pile-up control-distribution for data
		reco::VertexCollection vtx;
		reco::Vertex goodPV;
		unsigned int nGoodPV(0);
		fwlite::Handle<reco::VertexCollection> vtxHandle;
		vtxHandle.getByLabel(ev, "offlineSlimmedPrimaryVertices");
		if(vtxHandle.isValid() ) vtx = *vtxHandle;
		// Clean up vertex collection
		for(size_t ivtx=0; ivtx<vtx.size(); ++ivtx)
			{
			if(utils::isGoodVertex(vtx[ivtx]))
				{
				if(nGoodPV==0) goodPV=vtx[ivtx];
				nGoodPV++;
				}
			}

		// ----------------------------------------- Apply pileup reweighting
		// why don't use nGoodPV for Pile-Up?
		unsigned int num_inters = 0;
		if(isMC)
			{
			int ngenITpu = 0;
			fwlite::Handle < std::vector < PileupSummaryInfo > >puInfoH;
			puInfoH.getByLabel (ev, "slimmedAddPileupInfo");
			if (!puInfoH.isValid())
				{
				puInfoH.getByLabel( ev, "addPileupInfo" );
				if (!puInfoH.isValid()) {printf("collection PileupSummaryInfo with name slimmedAddPileupInfo or addPileupInfo does not exist\n"); exit(0);}
				}
			// so here we have valid puInfoH
			// otherwise exit was called
			for (std::vector < PileupSummaryInfo >::const_iterator it = puInfoH->begin (); it != puInfoH->end (); it++)
				{
				//if (it->getBunchCrossing () == 0) ngenITpu += it->getPU_NumInteractions ();
				// guys and Mara use getTrueNumInteractions :
				if (it->getBunchCrossing () == 0) ngenITpu += it->getTrueNumInteractions();
				}

			//ngenITpu = nGoodPV; // based on nvtx
			//puWeight = LumiWeights->weight (ngenITpu) * PUNorm[0];
			// So, in Pietro's approach ngenITpu is number of vertices in the beam crossing?
			//puWeight = direct_pileup_reweight[ngenITpu];
			// Mara does:
			//num_inters = puInfoH->at(0).getTrueNumInteractions(); // in 76 it seems to not work, returns 0 always
			// Using Pietro's PU number vertices:
			num_inters = ngenITpu;
			if (num_inters<100) {puWeight = direct_pileup_reweight[num_inters];}
			else {puWeight = 1.5e-16;}
			weight *= puWeight;
			// TODO: implement error margins of pile-up
			}
		else
			{
			// get data pile-up into num_inters
			// use number of good vertices for the data
			num_inters = nGoodPV;
			}

		// --------------- here the weighting/shaping of MC should be done
		// --------------------- save distributions of weights
		sum_weights += weight;
		sum_weights_raw += rawWeight;

		//num_inters = 1;
		if (num_inters>99) num_inters = 99;
		if (weightGen<0)
			{
			negative_event_nvtx[num_inters] += 1;
			negative_event_pernvtx_weight[num_inters] += weight;
			}
		else
			{
			positive_event_nvtx[num_inters] += 1;
			positive_event_pernvtx_weight[num_inters] += weight;
			}

		// -------------------------------------   Basic event selection

		// ---------------------- Orthogonalize Run2015B PromptReco+17Jul15 mix
		// let's remove Run2015B
		// if(isRun2015B)
		// {
		// if(!patUtils::exclusiveDataEventFilter(ev.eventAuxiliary().run(), isMC, isPromptReco ) ) continue;
		// }
		// it's not needed with the latest versions of RunB rereconstruction
		
		// -------------------------------------------------- Skip bad lumi
		// people say the new datasets for CMSSW76 don't have it implemented yet
		// testing if the procedure from 74 works with 76:
		if(!goodLumiFilter.isGoodLumi(ev.eventAuxiliary().run(),ev.eventAuxiliary().luminosityBlock())) continue; 
		
		// --------------------------------------------- apply trigger
		// ---------------- and require compatibilitiy of the event with the PD
		edm::TriggerResultsByName tr = ev.triggerResultsByName ("HLT");
		if (!tr.isValid ()){
			cout << "Trigger is not valid" << endl;
			return false;
			}

		if(debug){
			cout << "Printing trigger list" << endl;
			for(edm::TriggerNames::Strings::const_iterator trnames = tr.triggerNames().begin(); trnames!=tr.triggerNames().end(); ++trnames)
			cout << *trnames << endl;
			cout << "----------- End of trigger list ----------" << endl;
			//return 0;
		}

		// Need either to simulate the HLT (https://twiki.cern.ch/twiki/bin/view/CMS/TopTrigger#How_to_easily_emulate_HLT_paths) to match triggers.
		// Mara's triggers: HLT_Ele23_WPLoose_Gsf for electrons
		//                  HLT_IsoMu20 or HLT_IsoTkMu20 for muons
		bool eTrigger    ( utils::passTriggerPatterns(tr, "HLT_Ele23_WPLoose_Gsf*") );
		bool muTrigger   (
			utils::passTriggerPatterns (tr, "HLT_IsoMu20_v*", "HLT_IsoTkMu20_v*")
			);
		
		if(filterOnlySINGLEMU) {                    eTrigger = false; }
		if(filterOnlySINGLEE)  { muTrigger = false;                   }
		
		if (!(eTrigger || muTrigger)) continue;   //ONLY RUN ON THE EVENTS THAT PASS OUR TRIGGERS

		if(debug)
			{
			cout << "Set triggers" << endl;
			}

		// ------------------------------------------------- Apply MET filters
		//if( !isMC && !metFiler.passMetFilter( ev, isPromptReco)) continue;
		

		if(debug)
			{
			cout << "met filters are commented out here" << endl;
			}


		// ------------------------- event physics and the corresponding selection

		//------------------------- load all the objects we will need to access

		double rho = 0;
		fwlite::Handle<double> rhoHandle;
		rhoHandle.getByLabel(ev, "fixedGridRhoFastjetAll");
		if(rhoHandle.isValid() ) rho = *rhoHandle;
		
		reco::GenParticleCollection gen;
		fwlite::Handle<reco::GenParticleCollection> genHandle;
		genHandle.getByLabel(ev, "prunedGenParticles");
		if(genHandle.isValid() ) gen = *genHandle;
		
		// FIXME: Save time and don't load the rest of the objects when selecting by mctruthmode :)
		bool hasTop(false);
		int
		ngenLeptonsStatus3(0),
		ngenLeptonsNonTauSonsStatus3(0),
		ngenTausStatus3(0),
		ngenQuarksStatus3(0);
		//double tPt(0.), tbarPt(0.); // top pt reweighting - dummy value results in weight equal to 1 if not set in loop
		//float wgtTopPt(1.0), wgtTopPtUp(1.0), wgtTopPtDown(1.0);
		// TODO: what is this?? it should be the ttbar channel split
		// there was some vague answer from Pietro.....
		/*
		if(isMC)
			{
			// FIXME: Considering add support for different generators (based on PYTHIA6) for comparison.
			for(size_t igen=0; igen<gen.size(); igen++)
				{
				// FIXME: Should pass to the new status scheme from: https://github.com/cms-sw/cmssw/pull/7791
				// ////// if(!gen[igen].isHardProcess() && !gen[igen].isPromptFinalState()) continue;

				if(gen[igen].status() != 1 &&  gen[igen].status() !=2 && gen[igen].status() !=62 ) continue;
				int absid=abs(gen[igen].pdgId());
				// OK, so taus should be checked as status 2, and quarks as 71 or 23. More testing needed
				//if( absid==15 && hasWasMother(gen[igen]) ) cout << "Event " << iev << ", Particle " << igen << " has " << gen[igen].numberOfDaughters() << " daughters, pdgId " << gen[igen].pdgId() << " and status " << gen[igen].status() << ", mothers " << gen[igen].numberOfMothers() << ", pt " << gen[igen].pt() << ", eta " << gen[igen].eta() << ", phi " << gen[igen].phi() << ". isHardProcess is " << gen[igen].isHardProcess() << ", and isPromptFinalState is " << gen[igen].isPromptFinalState() << endl;


				////// if(absid==6 && gen[igen].isHardProcess()){ // particles of the hardest subprocess 22 : intermediate (intended to have preserved mass)
				if(absid==6 && gen[igen].status()==62)
					{
					// particles of the hardest subprocess 22 : intermediate (intended to have preserved mass). Josh says 62 (last in chain)
					hasTop=true;
					// FIXME: Top pT reweighting. 13 TeV values not propagated yet, so not using.
					//if(isTTbarMC){
					//  if(gen[igen].get("id") > 0) tPt=gen[igen].pt();
					//  else                        tbarPt=gen[igen].pt();
					//}
					} 


				//if(!gen[igen].isPromptFinalState() ) continue;
				if( (gen[igen].status() != 1 && gen[igen].status()!= 2 ) || !hasWasMother(gen[igen])) continue;

				if((absid==11 || absid==13) && hasLeptonAsDaughter(gen[igen]))
					cout << "Electron or muon " << igen << " has " << gen[igen].numberOfDaughters() << " daughter which is a lepton." << endl;

				if((absid==11 || absid==13) && gen[igen].status()==1)
					{
					ngenLeptonsStatus3++;

					if(!hasTauAsMother(gen[igen]))
						ngenLeptonsNonTauSonsStatus3++;
					}

				if(absid==15 && gen[igen].status()==2 )
					{
					ngenTausStatus3++; // This should be summed to ngenLeptonsStatus3 for the dilepton final states, not summed for the single lepton final states.
					//if(hasLeptonAsDaughter(gen[igen]))
					//	cout << "Tau " << igen << " has " << gen[igen].numberOfDaughters() << " daughter which is a lepton." << endl;
					}
				if (absid<=5) ngenQuarksStatus3++;
				}

			if(debug && (ngenTausStatus3==1 && ngenLeptonsStatus3==1 )  )
				cout << "Event: " << iev << ". Leptons: " << ngenLeptonsStatus3 << ". Leptons notaus: " << ngenLeptonsNonTauSonsStatus3 << ". Taus: " << ngenTausStatus3 << ". Quarks: " << ngenQuarksStatus3 << endl;

			// Dileptons:
			//    ttbar dileptons --> 1
			//    ttbar other     --> 2
			if(mctruthmode==1 && (ngenLeptonsStatus3+ngenTausStatus3!=2 || !hasTop )) continue;
			if(mctruthmode==2 && (ngenLeptonsStatus3+ngenTausStatus3==2 || !hasTop )) continue;
			// FIXME: port tt+bb splitting from 8 TeV (check the reference to the matched genjet)
			//if(mcTruthMode==1 && (ngenLeptonsStatus3!=2 || !hasTop || ngenBQuarksStatus23>=4)) continue;
			//if(mcTruthMode==2 && (ngenLeptonsStatus3==2 || !hasTop || ngenBQuarksStatus23>=4)) continue;
			//if(mcTruthMode==3 && (ngenBQuarksStatus23<4 || !hasTop))                           continue;

			// lepton-tau:
			//    ttbar ltau      --> 3
			//    ttbar dileptons --> 4
			//    ttbar ljets     --> 5
			//    ttbar hadrons   --> 6
			if(mctruthmode==3 && (ngenLeptonsNonTauSonsStatus3!=1 || ngenTausStatus3!=1  || !hasTop )) continue; // This is bugged, as it is obvious
			if(mctruthmode==4 && (ngenLeptonsNonTauSonsStatus3!=2                        || !hasTop )) continue;
			if(mctruthmode==5 && (ngenLeptonsNonTauSonsStatus3+ngenTausStatus3!=1        || !hasTop )) continue;

			bool isHad(false);
			if (
				(ngenLeptonsNonTauSonsStatus3!=1 || ngenTausStatus3!=1 ) &&
				(ngenLeptonsNonTauSonsStatus3!=2                      ) &&
				(ngenLeptonsNonTauSonsStatus3+ngenTausStatus3!=1      ) 
				)
			isHad=true;

			//if(mctruthmode==6 && (ngenLeptonsNonTauSonsStatus3!=0 || ngenTausStatus3!=0  || !hasTop )) continue;
			if(mctruthmode==6 && (!isHad || !hasTop )) continue;

			}

		if(debug) cout << "DEBUG: Event was not stopped by the ttbar sample categorization (either success, or it was not ttbar)" << endl;
		*/

		// FIXME: Top pT reweighting to be reactivated as soon as corrections are released
		// if(tPt>0 && tbarPt>0 && topPtWgt)
		//   {
		//   topPtWgt->computeWeight(tPt,tbarPt);
		//   topPtWgt->getEventWeight(wgtTopPt, wgtTopPtUp, wgtTopPtDown);
		//   wgtTopPtUp /= wgtTopPt;
		//   wgtTopPtDown /= wgtTopPt;
		//   }


		// ------------------------------------ actual particles?

		pat::MuonCollection muons;
		fwlite::Handle<pat::MuonCollection> muonsHandle;
		muonsHandle.getByLabel(ev, "slimmedMuons");
		if(muonsHandle.isValid() ) muons = *muonsHandle;

		pat::ElectronCollection electrons;
		fwlite::Handle<pat::ElectronCollection> electronsHandle;
		electronsHandle.getByLabel(ev, "slimmedElectrons");
		if(electronsHandle.isValid() ) electrons = *electronsHandle;

		// these "Collections" should be defined as:
		//typedef std::vector< pat::Electron > 	PatElectronCollection
 		// "define a PatElectronCollection as a vector of PatElectrons"


		/* No Jets, Photons, Taus and MET in Z
		pat::JetCollection jets;
		fwlite::Handle<pat::JetCollection>jetsHandle;
		jetsHandle.getByLabel(ev, "slimmedJets");
		if(jetsHandle.isValid() ) jets = *jetsHandle;

		pat::PhotonCollection photons;
		fwlite::Handle<pat::PhotonCollection> photonsHandle;
		photonsHandle.getByLabel(ev, "slimmedPhotons");
		if(photonsHandle.isValid() ) photons = *photonsHandle;

		pat::METCollection mets;
		fwlite::Handle<pat::METCollection> metsHandle;
		metsHandle.getByLabel(ev, "slimmedMETs");
		if(metsHandle.isValid() ) mets = *metsHandle;
		LorentzVector met = mets[0].p4 ();


		if(debug){
			// MET try:
			double mypt = mets[0].shiftedPt(pat::MET::METUncertainty::JetEnUp);
			cout << "MET = " << mets[0].pt() << ", JetEnUp: " << mypt << endl;
			LorentzVector myshiftedMet = mets[0].shiftedP4(pat::MET::METUncertainty::JetEnUp);
			cout << "MET = " << mets[0].pt() << ", JetEnUp: " << myshiftedMet.pt() << endl;
			}

		pat::TauCollection taus;
		fwlite::Handle<pat::TauCollection> tausHandle;
		tausHandle.getByLabel(ev, "slimmedTaus");
		if(tausHandle.isValid() ) taus = *tausHandle;
		*/

		//
		//
		// BELOW FOLLOWS THE ANALYSIS OF THE MAIN SELECTION WITH N-1 PLOTS. Whatever that means
		//
		//
		


		if(debug){
			cout << "got objects from the event, starting the analysis" << endl;
			}

		//
		// LEPTON ANALYSIS
		//
		
		// ------------------------------------ merging electrons and muons
		// Don't merge them in Z analysis
		// std::vector<patUtils::GenericLepton> leptons;
		// for(size_t l=0; l<electrons.size(); ++l) leptons.push_back(patUtils::GenericLepton (electrons[l] ));
		// for(size_t l=0; l<muons.size(); ++l)     leptons.push_back(patUtils::GenericLepton (muons[l]     ));

		// let's leave sort for now
		std::sort(electrons.begin(), electrons.end(), utils::sort_CandidatesByPt);
		std::sort(muons.begin(), muons.end(), utils::sort_CandidatesByPt);


		// ---------------------------------- leptons selection
		LorentzVector muDiff(0., 0., 0., 0.);

		std::vector<pat::Electron> selElectrons;
		std::vector<pat::Muon> selMuons;

		unsigned int nVetoE(0), nVetoMu(0);
		for(size_t ilep=0; ilep<muons.size (); ++ilep)
			{
			//patUtils::GenericLepton& muon = muons[ilep];
			pat::Muon& muon = muons[ilep];

			bool 
				passKin(true),     passId(true),     passIso(true),
				passVetoKin(true), passVetoId(true), passVetoIso(true);
						
			int lid(muon.pdgId());
						
			//apply muon corrections
			if(muCor)
				{
				TLorentzVector p4(muon.px(), muon.py(), muon.pz(), muon.energy());
				muCor->applyPtCorrection(p4, lid < 0 ? -1 : 1);
				if(isMC) muCor->applyPtSmearing(p4, lid < 0 ? -1 : 1, false);
				muDiff -= muon.p4();
				muon.setP4(LorentzVector(p4.Px(), p4.Py(), p4.Pz(), p4.E()));
				muDiff += muon.p4();
				}

			//no need for charge info any longer
			// lid = abs(lid);
			// TString lepStr(lid == 13 ? "mu" : "e");
					
			// no need to mess with photon ID // //veto nearby photon (loose electrons are many times photons...)
			// no need to mess with photon ID // double minDRlg(9999.);
			// no need to mess with photon ID // for(size_t ipho=0; ipho<selPhotons.size(); ipho++)
			// no need to mess with photon ID //   minDRlg=TMath::Min(minDRlg,deltaR(leptons[ilep].p4(),selPhotons[ipho].p4()));
			// no need to mess with photon ID // if(minDRlg<0.1) continue;

			// ---------------------------- kinematics
			double leta(fabs(muon.eta()));

			// ---------------------- Main lepton kin
			if(muon.pt() < 30.)                      passKin = false;
			if(leta > 2.1)                                    passKin = false;
			// if(lid == 11 && (leta > 1.4442 && leta < 1.5660)) passKin = false; // Crack veto

			// ---------------------- Veto lepton kin
			if (muon.pt () < 20)                      passVetoKin = false;
			if (leta > 2.1)                                    passVetoKin = false;
			// if (lid == 11 && (leta > 1.4442 && leta < 1.5660)) passVetoKin = false; // Crack veto

			//Cut based identification

			//std::vector<pat::Electron> dummyShit; dummyShit.push_back(leptons[ilep].el);

			// ------------------------- lepton IDs
			// passId     = lid == 11 ? patUtils::passId(electronVidMainId, myEvent, lepton.el) : patUtils::passId(lepton.mu, goodPV, patUtils::llvvMuonId::StdTight);
			// passVetoId = lid == 11 ? patUtils::passId(electronVidVetoId, myEvent, lepton.el) : patUtils::passId(lepton.mu, goodPV, patUtils::llvvMuonId::StdLoose);
			passId     = patUtils::passId(muon, goodPV, patUtils::llvvMuonId::StdTight);
			passVetoId = patUtils::passId(muon, goodPV, patUtils::llvvMuonId::StdLoose);
			// Notice:
			// any ID (Loose/Tight) means the lepton comes from "primary vertex"
			// (which is the vertex with max sum of Pt/Energy ot whatever else definition)
			// thus all leptons having some ID come from _the_same_ vertex

			// ------------------------- lepton isolation
			// passIso     = lid == 11 ? true : patUtils::passIso(lepton.mu, patUtils::llvvMuonIso::Tight); // Electron iso is included within the ID
			// passVetoIso = lid == 11 ? true : patUtils::passIso(lepton.mu, patUtils::llvvMuonIso::Loose); // Electron iso is included within the ID
			passIso     = patUtils::passIso(muon, patUtils::llvvMuonIso::Tight);
			passVetoIso = patUtils::passIso(muon, patUtils::llvvMuonIso::Loose);

			if     (passKin     && passId     && passIso)     selMuons.push_back(muon);
			else if(passVetoKin && passVetoId && passVetoIso) nVetoMu++;
			}


		for(size_t ilep=0; ilep<electrons.size (); ++ilep)
			{
			// patUtils::GenericLepton& muon = electrons[ilep];
			pat::Electron& electron = electrons[ilep];

			bool 
				passKin(true),     passId(true),     passIso(true),
				passVetoKin(true), passVetoId(true), passVetoIso(true);

			// ---------------------------- kinematics
			double leta(fabs(electron.superCluster()->eta()));

			// ---------------------- Main lepton kin
			if(electron.pt() < 30.)                      passKin = false;
			if(leta > 2.1)                                    passKin = false;
			if(leta > 1.4442 && leta < 1.5660) passKin = false; // Crack veto

			// ---------------------- Veto lepton kin
			if (electron.pt () < 20)                      passVetoKin = false;
			if (leta > 2.1)                                    passVetoKin = false;
			if (leta > 1.4442 && leta < 1.5660) passVetoKin = false; // Crack veto

			//Cut based identification

			//std::vector<pat::Electron> dummyShit; dummyShit.push_back(leptons[ilep].el);

			// ------------------------- lepton IDs
			// passId     = lid == 11 ? patUtils::passId(electronVidMainId, myEvent, lepton.el) : patUtils::passId(lepton.mu, goodPV, patUtils::llvvMuonId::StdTight);
			// passVetoId = lid == 11 ? patUtils::passId(electronVidVetoId, myEvent, lepton.el) : patUtils::passId(lepton.mu, goodPV, patUtils::llvvMuonId::StdLoose);
			passId     = patUtils::passId(electronVidMainId, myEvent, electron);
			passVetoId = patUtils::passId(electronVidVetoId, myEvent, electron);
			// Notice:
			// any ID (Loose/Tight) means the lepton comes from "primary vertex"
			// (which is the vertex with max sum of Pt/Energy ot whatever else definition)
			// thus all leptons having some ID come from _the_same_ vertex

			// ------------------------- lepton isolation
			// passIso     = lid == 11 ? true : patUtils::passIso(lepton.mu, patUtils::llvvMuonIso::Tight); // Electron iso is included within the ID
			// passVetoIso = lid == 11 ? true : patUtils::passIso(lepton.mu, patUtils::llvvMuonIso::Loose); // Electron iso is included within the ID
			passIso     = true;
			passVetoIso = true;

			if     (passKin     && passId     && passIso)     selElectrons.push_back(electron);
			else if(passVetoKin && passVetoId && passVetoIso) nVetoE++;
			}

		std::sort(selMuons.begin(),   selMuons.end(),   utils::sort_CandidatesByPt);
		std::sort(selElectrons.begin(),   selElectrons.end(),   utils::sort_CandidatesByPt);

		if (selMuons.size()<10) muon_selection_control[selMuons.size()] +=1;
		else muon_selection_control[10] +=1;

		if (selElectrons.size()<10) electron_selection_control[selElectrons.size()] +=1;
		else electron_selection_control[10] +=1;

		// LorentzVector recoMET = met;// FIXME REACTIVATE IT - muDiff;

		//
		// --------------------------------------------------
		//


		// find candidates among selected electrons and muons
		// go through the array of particles
		// select pairs which have the mass around Z mass
		std::vector<pair<int,int>> muon_z_pairs;
		if(selMuons.size()>=2)
			{
			for (int i = 0; i < selMuons.size() - 1; ++i)
				{
				int muon_i_id = selMuons[i].pdgId();
				for (int u = i + 1; u < selMuons.size(); ++u)
					{
					int muon_u_id = selMuons[u].pdgId();
					Double_t dileptonSystem_mass = (selMuons[u].p4() + selMuons[i].p4()).M();
					// different sign and mass in the range
					if ( (muon_u_id*muon_i_id < 0) && dileptonSystem_mass > 50  && dileptonSystem_mass < 150 )
						{
						std::pair<int,int> good_pair;
						good_pair.first =  i;
						good_pair.second = u;
						muon_z_pairs.push_back(good_pair);
						}
					}
				}
			}

		n_good_muon_pairs += muon_z_pairs.size();

		std::vector<pair<int,int>> electron_z_pairs;
		if (selElectrons.size()>=2)
			{
			for (int i = 0; i < selElectrons.size() - 1; ++i)
				{
				int i_id = selElectrons[i].pdgId();
				for (int u = i + 1; u < selElectrons.size(); ++u)
					{
					int u_id = selElectrons[u].pdgId();
					Double_t dileptonSystem_mass = (selElectrons[u].p4() + selElectrons[i].p4()).M();
					// different sign and mass in the range
					if ( (u_id*i_id < 0) && dileptonSystem_mass > 50.  && dileptonSystem_mass < 150. )
						{
						std::pair<int,int> good_pair;
						good_pair.first =  i;
						good_pair.second = u;
						electron_z_pairs.push_back(good_pair);
						}
					}
				}
			}

		n_good_electron_pairs += electron_z_pairs.size();

		if(debug){
			cout << "Finished processing first event in the first file, exiting" << endl;
			//return 0;
			break;
			}

		} // End single file event loop

	fprintf(csv_out, "N good muon pairs:%d\n",     n_good_muon_pairs);
	fprintf(csv_out, "distr of N selected muons:%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n", muon_selection_control[0],muon_selection_control[1],muon_selection_control[2],muon_selection_control[3],muon_selection_control[4], muon_selection_control[5],muon_selection_control[6],muon_selection_control[7],muon_selection_control[8],muon_selection_control[9], muon_selection_control[10]);
	fprintf(csv_out, "N good electron pairs:%d\n", n_good_electron_pairs);
	fprintf(csv_out, "distr of N selected electrons:%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n", electron_selection_control[0],electron_selection_control[1],electron_selection_control[2],electron_selection_control[3],electron_selection_control[4], electron_selection_control[5],electron_selection_control[6],electron_selection_control[7],electron_selection_control[8],electron_selection_control[9], electron_selection_control[10]);

	printf("\n");

	delete file;
	} // End loop on files

fclose(csv_out);

// if(saveSummaryTree)
// 	{
// 	TDirectory* cwd = gDirectory;
// 	summaryFile->cd();
// 	summaryTree->Write();
// 	summaryFile->Close();
// 	delete summaryFile;
// 	cwd->cd();
// 	}


if(nMultiChannel>0) cout << "Warning! There were " << nMultiChannel << " multi-channel events out of " << totalEntries << " events!" << endl;
printf ("\n");

//##############################################
//########     SAVING HISTO TO FILE     ########
//##############################################
//save control plots to file
printf ("Results save in %s\n", outUrl.Data());

//save all to the file
// TFile *ofile = TFile::Open (outUrl + ".root", "recreate");


// ofile->Close();

if (outTxtFile) fclose (outTxtFile);

// Now that everything is done, dump the list of lumiBlock that we processed in this job
if(!isMC){
	// FIXME: when lumi certificate is ready for rereco data, check that these work
	goodLumiFilter.FindLumiInFiles(urls);
	goodLumiFilter.DumpToJson(((outUrl.ReplaceAll(".root",""))+".json").Data());
	}

}

