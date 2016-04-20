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
		//
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
					for(size_t ilep=0; ilep<leptons.size(); ilep++) {
						LorentzVector lepton = leptons[ilep].p4();
						double varSign( (ivar==LESUP ? 1.0 : (ivar==LESDOWN ? -1.0 : 0.0) ) );
						int id( abs(leptons[ilep].pdgId()) );
						double sf(1.0);
						if(id==13) sf=(1.0+varSign*0.01);
						if(id==11) {
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
bool passPFJetID(std::string label,
								 pat::Jet jet){
	
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
	if(label=="Loose")
		passID = ( ((nhf<0.99 && nef<0.99 && nconst>1) && ((abs(eta)<=2.4 && chf>0 && nch>0 && cef<0.99) || abs(eta)>2.4)) && abs(eta)<=3.0 );
	if(label=="Tight")
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
	while ((part->numberOfDaughters()>0)) {
		const reco::Candidate* DaughterPart = part->daughter(0);
		// cout << "\t\t Daughter: " << DaughterPart->pdgId() << " with status " << DaughterPart->status() << endl;
		if (fabs(DaughterPart->pdgId()) == 11 || fabs(DaughterPart->pdgId() == 13)){
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
	while ((part->numberOfMothers()>0)) {
		const reco::Candidate* MomPart =part->mother();
		if (fabs(MomPart->pdgId())==24){
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
	if(p.numberOfMothers()==0) return foundTau;
	const reco::Candidate* part = &p; //(p.mother());
	// loop on the mother particles to check if it has a tau as mother
	while ((part->numberOfMothers()>0)) {
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
			//      systVars.push_back(); systVars.push_back();

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
	//    KEY: TTreeMetaData;1
	//                         KEY: TTreeParameterSets;1
	//                                                   KEY: TTreeParentage;1
	//                                                                         KEY: TTreeEvents;1
	//                                                                                            KEY: TTreeLuminosityBlocks;1
	//                                                                                                                         KEY: TTreeRuns;
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
	MuScleFitCorrector *muCor = NULL; // FIXME: MuScle fit corrections for 13 TeV not available yet (more Zs are needed) getMuonCorrector (jecDir, dtag);

	// --------------------------------------- lepton efficiencies
	LeptonEfficiencySF lepEff;
	
	// --------------------------------------- b-tagging 
	// Prescriptions taken from: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation74X

	// b-tagging working points for 50ns 
	//   (pfC|c)ombinedInclusiveSecondaryVertexV2BJetTags
	//      v2CSVv2L 0.605
	//      v2CSVv2M 0.890
	//      v2CSVv2T 0.970
	double
		btagLoose(0.605),
		btagMedium(0.890),
		btagTight(0.970);

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


	TString
		electronIdMainTag("cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
		electronIdVetoTag("cutBasedElectronID-Spring15-25ns-V1-standalone-tight");

	// --------------------------------------- pileup weighting
	/*
	edm::LumiReWeighting * LumiWeights = NULL;
	utils::cmssw::PuShifter_t PuShifters;
	double PUNorm[] = { 1, 1, 1 };
	if (isMC)
		{
			std::vector<double> dataPileupDistributionDouble = runProcess.getParameter < std::vector < double >>("datapileup");
			std::vector<float> dataPileupDistribution;
			for (unsigned int i = 0; i < dataPileupDistributionDouble.size (); i++)
				{
					dataPileupDistribution.push_back (dataPileupDistributionDouble[i]);
				}
			std::vector<float> mcPileupDistribution;
			utils::getMCPileupDistributionFromMiniAOD(urls, dataPileupDistribution.size (), mcPileupDistribution);
			while(mcPileupDistribution.size() < dataPileupDistribution.size()) mcPileupDistribution.push_back(0.0);
			while(mcPileupDistribution.size() > dataPileupDistribution.size()) dataPileupDistribution.push_back(0.0);
			gROOT->cd ();             //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE
			LumiWeights = new edm::LumiReWeighting(mcPileupDistribution, dataPileupDistribution);
			PuShifters = utils::cmssw::getPUshifters(dataPileupDistribution, 0.05);
			utils::getPileupNormalization(mcPileupDistribution, PUNorm, LumiWeights, PuShifters);
		}
	 */
	 // pile-up is done directly with direct_pileup_reweight
	
	gROOT->cd ();                 //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE
	
	//higgs::utils::EventCategory eventCategoryInst(higgs::utils::EventCategory::EXCLUSIVE2JETSVBF); //jet(0,>=1)+vbf binning

	std::vector<double> direct_pileup_reweight = runProcess.getParameter < std::vector < double >>("pileup_reweight_direct");
 

	// --------------------------------------- hardcoded MET filter
	patUtils::MetFilter metFiler;
	if(!isMC) { 
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

	// removing the SmartSelectionMonitor
	// SmartSelectionMonitor mon;


	TH1D* singlelep_ttbar_initialevents  = (TH1D*) new TH1D("singlelep_ttbar_init",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  ); 
	TH1D* singlelep_ttbar_preselectedevents = (TH1D*) new TH1D("singlelep_ttbar_presele",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  ); 
	TH1D* singlelep_ttbar_selected_mu_events = (TH1D*) new TH1D("singlelep_ttbar_sele_mu",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  ); 
	TH1D* singlelep_ttbar_selected_el_events = (TH1D*) new TH1D("singlelep_ttbar_sele_el",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  ); 
	TH1D* singlelep_ttbar_selected2_mu_events = (TH1D*) new TH1D("singlelep_ttbar_sele2_mu",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  ); 
	TH1D* singlelep_ttbar_selected2_el_events = (TH1D*) new TH1D("singlelep_ttbar_sele2_el",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  ); 

	TH1D* singlelep_ttbar_maraselected_mu_events = (TH1D*) new TH1D("singlelep_ttbar_marasele_mu",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  ); 
	TH1D* singlelep_ttbar_maraselected_el_events = (TH1D*) new TH1D("singlelep_ttbar_marasele_el",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  ); 

/*
	// ensure proper normalization
	TH1D* normhist = (TH1D*) mon.addHistogram(new TH1D("initNorm", ";;Events", 5, 0., 5.));
	normhist->GetXaxis()->SetBinLabel (1, "Gen. Events");
	normhist->GetXaxis()->SetBinLabel (2, "Events");
	normhist->GetXaxis()->SetBinLabel (3, "PU central");
	normhist->GetXaxis()->SetBinLabel (4, "PU up");
	normhist->GetXaxis()->SetBinLabel (5, "PU down");

	//event selection - charged Higgs
	TH1D* h = (TH1D*) mon.addHistogram (new TH1D ("chhiggseventflowdilep", ";;Events", 6, 0., 6.));
	h->GetXaxis()->SetBinLabel (1, "#geq 2 iso leptons");
	h->GetXaxis()->SetBinLabel (2, "M_{ll} veto");
	h->GetXaxis()->SetBinLabel (3, "#geq 2 jets");
	h->GetXaxis()->SetBinLabel (4, "E_{T}^{miss}");
	h->GetXaxis()->SetBinLabel (5, "op. sign");
	h->GetXaxis()->SetBinLabel (6, "#geq 2 b-tags");
	h = (TH1D*) mon.addHistogram(new TH1D ("chhiggsalteventflowdilep", ";;Events", 8, 0., 8.));
	h->GetXaxis()->SetBinLabel (1, "#geq 2 iso leptons");
	h->GetXaxis()->SetBinLabel (2, "M_{ll} veto");
	h->GetXaxis()->SetBinLabel (3, "#geq 2 jets");
	h->GetXaxis()->SetBinLabel (4, "E_{T}^{miss}");
	h->GetXaxis()->SetBinLabel (5, "op. sign");
	h->GetXaxis()->SetBinLabel (6, "= 0 b-tags");
	h->GetXaxis()->SetBinLabel (7, "= 1 b-tags");
	h->GetXaxis()->SetBinLabel (8, "#geq 2 b-tags");
	h = (TH1D*) mon.addHistogram(new TH1D ("chhiggseventflowslep", ";;Events", 6, 0., 6.));
	h->GetXaxis()->SetBinLabel (1, "1 iso lepton");
	h->GetXaxis()->SetBinLabel (2, "#geq 2 jets");
	h->GetXaxis()->SetBinLabel (3, "E_{T}^{miss}");
	h->GetXaxis()->SetBinLabel (4, "#geq 1 b-tag");
	h->GetXaxis()->SetBinLabel (5, "1 #tau");
	h->GetXaxis()->SetBinLabel (6, "op. sign");
	h = (TH1D*) mon.addHistogram(new TH1D("chhiggsalteventflowslep", ";;Events", 7, 0., 7.));
	h->GetXaxis()->SetBinLabel (1, "1 iso lepton");
	h->GetXaxis()->SetBinLabel (2, "E_{T}^{miss}");
	h->GetXaxis()->SetBinLabel (3, "1 #tau");
	h->GetXaxis()->SetBinLabel (4, "op. sign");
	h->GetXaxis()->SetBinLabel (5, "= 0 b-tags");
	h->GetXaxis()->SetBinLabel (6, "= 1 b-tags");
	h->GetXaxis()->SetBinLabel (7, "#geq 2 b-tags");

	// event selection - cross section
	h = (TH1D*) mon.addHistogram (new TH1D ("xseceventflowdilep", ";;Events", 6, 0., 6.));
	h->GetXaxis()->SetBinLabel (1, "#geq 2 iso leptons");
	h->GetXaxis()->SetBinLabel (2, "M_{ll} veto");
	h->GetXaxis()->SetBinLabel (3, "#geq 2 jets");
	h->GetXaxis()->SetBinLabel (4, "E_{T}^{miss}");
	h->GetXaxis()->SetBinLabel (5, "op. sign");
	h->GetXaxis()->SetBinLabel (6, "#geq 2 b-tags");
	h = (TH1D*) mon.addHistogram(new TH1D ("xsecalteventflowdilep", ";;Events", 8, 0., 8.));
	h->GetXaxis()->SetBinLabel (1, "#geq 2 iso leptons");
	h->GetXaxis()->SetBinLabel (2, "M_{ll} veto");
	h->GetXaxis()->SetBinLabel (3, "#geq 2 jets");
	h->GetXaxis()->SetBinLabel (4, "E_{T}^{miss}");
	h->GetXaxis()->SetBinLabel (5, "op. sign");
	h->GetXaxis()->SetBinLabel (6, "= 0 b-tags");
	h->GetXaxis()->SetBinLabel (7, "= 1 b-tags");
	h->GetXaxis()->SetBinLabel (8, "#geq 2 b-tags");
	h = (TH1D*) mon.addHistogram(new TH1D ("xseceventflowslep", ";;Events", 6, 0., 6.));
	h->GetXaxis()->SetBinLabel (1, "1 iso lepton");
	h->GetXaxis()->SetBinLabel (2, "#geq 2 jets");
	h->GetXaxis()->SetBinLabel (3, "E_{T}^{miss}");
	h->GetXaxis()->SetBinLabel (4, "#geq 1 b-tag");
	h->GetXaxis()->SetBinLabel (5, "1 #tau");
	h->GetXaxis()->SetBinLabel (6, "op. sign");
	h = (TH1D*) mon.addHistogram(new TH1D("xsecalteventflowslep", ";;Events", 7, 0., 7.));
	h->GetXaxis()->SetBinLabel (1, "1 iso lepton");
	h->GetXaxis()->SetBinLabel (2, "E_{T}^{miss}");
	h->GetXaxis()->SetBinLabel (3, "1 #tau");
	h->GetXaxis()->SetBinLabel (4, "op. sign");
	h->GetXaxis()->SetBinLabel (5, "= 0 b-tags");
	h->GetXaxis()->SetBinLabel (6, "= 1 b-tags");
	h->GetXaxis()->SetBinLabel (7, "#geq 2 b-tags");

	h = (TH1D*) mon.addHistogram( new TH1D("nvtx_pileup"         , ";;Events", 100, 0., 100.));
	h = (TH1D*) mon.addHistogram( new TH1D("nvtx_singlemu_pileup", ";;Events", 100, 0., 100.));
	h = (TH1D*) mon.addHistogram( new TH1D("nvtx_singlee_pileup" , ";;Events", 100, 0., 100.));

	
	// Setting up control categories to be analyzed
	std::vector < TString > controlCats;
	controlCats.clear ();
	controlCats.push_back("step1");
	controlCats.push_back("step2");
	controlCats.push_back("step3");
	controlCats.push_back("step4");
	controlCats.push_back("step5");
	controlCats.push_back("step6");
	
	controlCats.push_back("altstep1");
	controlCats.push_back("altstep2");
	controlCats.push_back("altstep3");
	controlCats.push_back("altstep4");
	controlCats.push_back("altstep5");
	controlCats.push_back("altstep6");
	controlCats.push_back("altstep7");
	controlCats.push_back("altstep8");

	for (size_t k = 0; k < controlCats.size (); ++k)
		{
			TString icat (controlCats[k]);

			//pu control to be completed
			mon.addHistogram (new TH1D (icat+"nvtx",    ";Vertices;Events", 100, 0., 100.));
			mon.addHistogram (new TH1D (icat+"nvtxraw", ";Vertices;Events", 100, 0., 100.));
			mon.addHistogram (new TH1D (icat+"rho", ";#rho;Events"        ,  50, 0.,  25.));
			

			//tau control to be completed
			TH1* htaus = mon.addHistogram (new TH1D (icat + "ntaus", ";Tau multiplicity;Events", 5, 0., 5.));
			for (int ibin = 1; ibin <= htaus->GetXaxis ()->GetNbins (); ibin++)
				{
					TString label ("");
					if (ibin == h->GetXaxis ()->GetNbins ())
						label += "#geq";
					else
						label += "=";
					label += (ibin - 1);
					htaus->GetXaxis ()->SetBinLabel (ibin, label);
				}
			mon.addHistogram( new TH1D(icat+"tauleadpt",              ";p_{T}^{#tau};Events", 30,  0.,  300.  ));
			mon.addHistogram( new TH1D(icat+"tauleadeta",             ";#eta^{#tau};Events",  50, -2.6,   2.6 ));
			mon.addHistogram( new TH1D(icat+"tauinclusivept",         ";p_{T}^{#tau};Events", 30,  0.,  300.  ));
			mon.addHistogram( new TH1D(icat+"tauinclusiveeta",        ";#eta^{#tau};Events",  50, -2.6,   2.6 ));
			mon.addHistogram( new TH1D(icat+"tauinclusivecharge",     ";p_{T}^{#tau};Events",  5, -2.,    2.  ));
			mon.addHistogram( new TH1D(icat+"tauinclusivedz",         ";dz^{#tau};Events",    50,  0.,   10.  ));
			mon.addHistogram( new TH1D(icat+"tauinclusivevz",         ";vz^{#tau};Events",    50,  0.,   10.  ));
			mon.addHistogram( new TH1D(icat+"tauinclusiveemfraction", ";emf^{#tau};Events",   50,  0.,    5.  ));
			mon.addHistogram( new TH1D(icat+"tauinclusivedizeta"    , ";dZ^{#tau};Events",    50,  0.,   10.  ));

			//lepton control
			mon.addHistogram( new TH1D(icat+"inclusivept",      ";Transverse momentum [GeV];Events",             50, 0.,  500.  ));
			mon.addHistogram( new TH1D(icat+"leadleptonpt",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  )); 
			mon.addHistogram( new TH1D(icat+"leadleptoneta",    ";Pseudo-rapidity;Events",                       50, 0.,    2.5 ));
			mon.addHistogram( new TH1D(icat+"trailerleptonpt",  ";Transverse momentum [GeV];Events",             50, 0.,  500.  ));
			mon.addHistogram( new TH1D(icat+"trailerleptoneta", ";Pseudo-rapidity;Events",                       50, 0.,    2.5 ));
			mon.addHistogram( new TH1D(icat+"pte",              ";Electron transverse momentum [GeV];Events",    50, 0.,  500.  ));
			mon.addHistogram( new TH1D(icat+"ptmu",             ";Muon transverse momentum [GeV];Events",        50, 0.,  500.  ));
			mon.addHistogram( new TH1D(icat+"qt",               ";Transverse momentum [GeV];Events / (1 GeV)", 1500, 0., 1500.  ));
			mon.addHistogram( new TH1D(icat+"emva", "; e-id MVA; Electrons", 50, 0.95,1.0) );
			
			// Dilepton control
			mon.addHistogram( new TH1D(icat+"sumptll",      ";Sum of lepton transverse momenta [GeV];Events",                     75, 0.,  750.  ));
			mon.addHistogram( new TH1D(icat+"mll",          ";Dilepton invariant mass [GeV];Events",                              50, 0.,  500.  ));
			mon.addHistogram( new TH1D(icat+"ptll",         ";Dilepton transverse momentum [GeV];Events",                         75, 0.,  750.  ));
			mon.addHistogram( new TH1D(icat+"yll",          ";Rapidity;Events",                                                   50, 0.,    3.  ));
			mon.addHistogram( new TH1D(icat+"dilarccosine", ";#theta(l,l') [rad];Events",                                         64, 0.,    3.2 ));

			mon.addHistogram( new TH1D(icat+"mtsum",        ";M_{T}(l^{1},E_{T}^{miss})+M_{T}(l^{2},E_{T}^{miss}) [GeV];Events", 100, 0.,   10.  ));
			mon.addHistogram( new TH1D(icat+"ht",           ";H_{T} [GeV];Events",                                                75, 0., 1500.  ));
			mon.addHistogram( new TH1D(icat+"htb",          ";H_{T} (bjets) [GeV];Events",                                        75, 0., 1500.  ));
			mon.addHistogram( new TH1D(icat+"htnol",        "; H_[T] (no leptons) [GeV];Events",                                  75, 0., 1500.  ));
			mon.addHistogram( new TH1D(icat+"htbnol",       "; H_[T] (bjets, no leptons) [GeV];Events",                           75, 0., 1500.  ));

			// add plots for extra leptons in the event
			// add plots for third lepton pt etc

			mon.addHistogram (new TH1D(icat+"csv",           ";Combined Secondary Vertex;Jets",    50, 0.,    1. ));
			mon.addHistogram (new TH1D(icat+"csvb",          ";Combined Secondary Vertex;Jets",    50, 0.,    1. ));
			mon.addHistogram (new TH1D(icat+"csvc",          ";Combined Secondary Vertex;Jets",    50, 0.,    1. ));
			mon.addHistogram (new TH1D(icat+"csvothers",     ";Combined Secondary Vertex;Jets",    50, 0.,    1. ));

			mon.addHistogram (new TH1D(icat+"leadjetpt",     ";Transverse momentum [GeV];Events", 100, 0., 1000. ));
			mon.addHistogram (new TH1D(icat+"trailerjetpt",  ";Transverse momentum [GeV];Events", 100, 0., 1000. ));
			mon.addHistogram (new TH1D(icat+"fwdjeteta",     ";Pseudo-rapidity;Events",            60, 0.,    3. ));
			mon.addHistogram (new TH1D(icat+"leadjeteta",    ";Pseudo-rapidity;Events",            60, 0.,    3. ));
			mon.addHistogram (new TH1D(icat+"trailerjeteta", ";Pseudo-rapidity;Events",            60, 0.,    3. ));
			mon.addHistogram (new TH1D(icat+"cenjeteta",     ";Pseudo-rapidity;Events",            60, 0.,    3. ));
			
			TH1* hbtags   = mon.addHistogram(new TH1D(icat+"nbtags",   ";b-tag multiplicity;Events", 5, 0., 5. ));
			TH1* hjets    = mon.addHistogram(new TH1D(icat+"njets",    ";Jet multiplicity;Events",   6, 0., 6. ));
			for (int ibin = 1; ibin <= hjets->GetXaxis ()->GetNbins (); ibin++)
				{
					TString label ("");
					if (ibin == h->GetXaxis()->GetNbins() || (TString(h->GetName()).Contains("btags") && ibin == h->GetXaxis()->GetNbins()-1 ) )
						label += "#geq";
					else
						label += "=";
					label += (ibin - 1);
					hjets   ->GetXaxis()->SetBinLabel(ibin, label);
					hbtags  ->GetXaxis()->SetBinLabel(ibin, label);
				}
		 
			mon.addHistogram (new TH1D (icat + "mindphijmet",    ";min #Delta#phi(jet,E_{T}^{miss});Events",     40,    0.,    4.  ));
			mon.addHistogram (new TH1D (icat + "mindphijmetNM1", ";min #Delta#phi(jet,E_{T}^{miss});Events",     40,    0.,    4.  ));
			mon.addHistogram (new TH1D (icat + "balance",        ";E_{T}^{miss}/q_{T};Events",                   25,    0.,    2.5 ));
			mon.addHistogram (new TH1D (icat + "balanceNM1",     ";E_{T}^{miss}/q_{T};Events",                   25,    0.,    2.5 ));
			mon.addHistogram (new TH1D (icat + "axialmet",       ";Axial missing transvere energy [GeV];Events", 50, -100.,  400.  ));
			mon.addHistogram (new TH1D (icat + "axialmetNM1",    ";Axial missing transvere energy [GeV];Events", 50, -100.,  400.  ));
			mon.addHistogram (new TH1D (icat + "met",            ";Missing transverse energy [GeV];Events",      50,    0., 1000.  ));
			mon.addHistogram (new TH1D (icat + "recomet",        ";Missing transverse energy [GeV];Events",      50,    0., 1000.  ));
			mon.addHistogram (new TH1D (icat + "mt",             ";Transverse mass;Events",                      50,    0.,  500.  ));
			mon.addHistogram (new TH1D (icat + "mtresponse",     ";Transverse mass response;Events",            100,    0.,    2.  ));
			mon.addHistogram (new TH1D (icat + "mtcheckpoint",   ";Transverse mass [GeV];Events",               160,  150., 1750.  ));
			mon.addHistogram (new TH1D (icat + "metcheckpoint",  ";Missing transverse energy [GeV];Events",     100,    0.,  500.  ));

		} // End of loop on controlCats


	//
	// STATISTICAL ANALYSIS
	//
	TH1D* Hoptim_systs = (TH1D*) mon.addHistogram (new TH1D ("optim_systs", ";syst;", nSystVars, 0, nSystVars));
	for (size_t ivar=0; ivar<nSystVars; ++ivar) Hoptim_systs->GetXaxis()->SetBinLabel(ivar+1, systVars[ivar]);

	// Final distributions to compute systematics on
	for(size_t ivar=0; ivar<nSystVars; ++ivar) 
		{
			TString var=systVars[ivar];
			
			// dilepton or both
			mon.addHistogram(new TH1D("finalnbjets"         +var, ";b-jet multiplicity;Events",     6, 0.,   6. ));        
			mon.addHistogram(new TH1D("finalmt"             +var, ";Transverse mass [GeV];Events", 50, 0., 500. ));
			
			// lepton-tau
			mon.addHistogram(new TH1D("finaltaur"           +var, ";R^{#tau};Events",                    10,  0.,   1.   ));
			mon.addHistogram(new TH1D("finaltaupolarization"+var, ";Y^{#tau};Events",                    40, -1.,   3.   ));
			mon.addHistogram(new TH1D("finaldphilepmet"     +var, ";#Delta#phi(#tau_{h}-#it{l});Events", 60,  0.  , 3.15 ));
			mon.addHistogram(new TH1D("finaldphitaumet"     +var, ";#Delta#phi(#tau_{h}-MET);Events",    60,  0.,   3.15 ));
			mon.addHistogram(new TH1D("finaldphileptau"     +var, ";#Delta#phi(#it{l}-#tau_{h});Events", 60,  0.,   3.15 ));
			mon.addHistogram(new TH1D("finaltaupt"          +var, ";p_{T}^{#tau} [GeV];Events",          50,  0., 500.   ));
			mon.addHistogram(new TH1D("finalmutaumass"      +var, ";M_{#mu#tau_{h}} [GeV];Events",       50,  0., 500.   ));

		}

*/

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
//csv_out = fopen(const char *filename, const char *mode);
string FileName = ((outUrl.ReplaceAll(".root",""))+".csv").Data();
//const char file_name = FileName.c_str();
csv_out = fopen(FileName.c_str(), "w");
fprintf(csv_out, "Headers\n");
//fprintf(csv_out, "kino:\npl.E, plb.E, pb.E, pbb.E,\nprest-sqr, prest-XY-sqr, met.pt,\nprest-o-plpb, pl-o-pb, plb-o-pbb,\nsame 3 angles\n");
//fprintf(csv_out, "kino2:\n(pl1) x,y,z,e\n(pl2) x,y,z,e\n(pb1) x,y,z,e\n(pb2) x,y,z,e\n\n");

fprintf(csv_out, "acceptances:filename, num_events, sum_rawWeight, sum_weight, cross_sum_rawWeight,cross_sum_weight, oursel_sum_rawWeight,oursel_sum_weight, marasel_sum_rawWeight,marasel_sum_weight\n");

fprintf(csv_out, "crossel:pu_num_inters,nGoodPV, rawWeight, weight, isElectron, l_px,l_py,l_pz,l_e, b1_px,b1_py,b1_pz,b1_e, j1_px,j1_py,j1_pz,j1_e,j2_px,j2_py,j2_pz,j2_e\n");
fprintf(csv_out, "oursel:pu_num_inters,nGoodPV, rawWeight, weight, isElectron, l_px,l_py,l_pz,l_e, tau_px,tau_py,tau_pz,tau_e, b1_px,b1_py,b1_pz,b1_e, j1_px,j1_py,j1_pz,j1_e,j2_px,j2_py,j2_pz,j2_e\n");
fprintf(csv_out, "marasel:pu_num_inters,nGoodPV, rawWeight, weight, isElectron, l_px,l_py,l_pz,l_e, b1_px,b1_py,b1_pz,b1_e,b2_px,b2_py,b2_pz,b2_e, j1_px,j1_py,j1_pz,j1_e,j2_px,j2_py,j2_pz,j2_e,j3_px,j3_py,j3_pz,j3_e,j4_px,j4_py,j4_pz,j4_e\n");

for(size_t f=0; f<urls.size();++f){
	fprintf(csv_out, "Processing file: %s\n", urls[f].c_str());
	TFile* file = TFile::Open(urls[f].c_str());
	fwlite::Event ev(file);
	printf ("Scanning the ntuple %2lu/%2lu :",f+1, urls.size());

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
	for (int i=0; i<100; i++)
		{
		negative_event_nvtx[i] = 0;
		positive_event_nvtx[i] = 0;
		}

	int treeStep (ev.size()/50);
	//DuplicatesChecker duplicatesChecker;
	//int nDuplicates(0);
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


		// ---------------------------------- these are weird NLO -1 weights
		// TODO: figure out what are these really
		// Take into account the negative weights from some NLO generators (otherwise some phase space will be double counted)
		double weightGen(1.);
		if(isNLOMC)
			{
			//double weightGen(0.);
			//double weightLhe(0.);

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

		// ---------------------------------- pileup weight
		double weight           (1.0);
		double rawWeight        (1.0);
		double TotalWeight_plus (1.0);
		double TotalWeight_minus(1.0);
		double puWeight         (1.0);


		// This must remain deactivated if you use HT-binned samples (it was for pthat-binned samples)
		// if (isV0JetsMC)
		//   {
		//     fwlite::Handle < LHEEventProduct > lheEPHandle;
		//     lheEPHandle.getByLabel (ev, "externalLHEProducer");
		//     mon.fillHisto ("nup", "", lheEPHandle->hepeup ().NUP, 1);
		//     if (lheEPHandle->hepeup ().NUP > 5)  continue;
		//     mon.fillHisto ("nupfilt", "", lheEPHandle->hepeup ().NUP, 1);
		//   }

		// HT-binned samples stitching: https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorking2015#MC_and_data_samples
		// -------------- it should be the creeppish merging of LO and NLO sets
		// not used now at all?
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
			for ( size_t idxParticle = 0; idxParticle < numParticles; ++idxParticle ) {
				int absPdgId = TMath::Abs(lheEvent.IDUP[idxParticle]);
				int status = lheEvent.ISTUP[idxParticle];
				if ( status == 1 && ((absPdgId >= 1 && absPdgId <= 6) || absPdgId == 21) ) {
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
			for (std::vector < PileupSummaryInfo >::const_iterator it = puInfoH->begin (); it != puInfoH->end (); it++)
				{
				if (it->getBunchCrossing () == 0) ngenITpu += it->getPU_NumInteractions ();
				}

			//ngenITpu = nGoodPV; // based on nvtx
			//puWeight = LumiWeights->weight (ngenITpu) * PUNorm[0];
			// So, in Pietro's approach ngenITpu is number of vertices in the beam crossing?
			//puWeight = direct_pileup_reweight[ngenITpu];
			// Mara does:
			//num_inters = puInfoH->at(0).getTrueNumInteractions(); // in 76 it seems to not work, returns 0 always
			// Using Pietro's PU number vertices:
			num_inters = ngenITpu;
			// FIXME: hopefully the length of the array is enough. Increse to 100 bins.
			if (num_inters<40) {puWeight = direct_pileup_reweight[num_inters];}
			else {puWeight = 1.5e-16;}
			weight *= puWeight;//Weight; //* puWeight;
			// implement error margins of pile-up
			//TotalWeight_plus =  PuShifters[utils::cmssw::PUUP]  ->Eval (ngenITpu) * (PUNorm[2]/PUNorm[0]);
			//TotalWeight_minus = PuShifters[utils::cmssw::PUDOWN]->Eval (ngenITpu) * (PUNorm[1]/PUNorm[0]);
			}

		// --------------- here the weighting/shaping of MC should be done
		// --------------------- save distributions of weights
		sum_weights += weight;
		sum_weights_raw += rawWeight;

		if (weightGen<0) negative_event_nvtx[num_inters] += 1;
                else positive_event_nvtx[num_inters] += 1;

		// TODO: make separate weight and number of events distrobutions
		//mon.fillHisto("initNorm", tags, 0., weightGen); // Should be all 1, but for NNLO samples there are events weighting -1
		//mon.fillHisto("initNorm", tags, 1., weightGen); // Should be all 1, but for NNLO samples there are events weighting -1
		//mon.fillHisto("initNorm", tags, 2., puWeight);
		//mon.fillHisto("initNorm", tags, 3., TotalWeight_plus);
		//mon.fillHisto("initNorm", tags, 4., TotalWeight_minus);
		// probably, these are N events after re-forming MC
		
		// ############################################   EVENT LOOP STARTS

		// ---------------------- Orthogonalize Run2015B PromptReco+17Jul15 mix
		// let's remove Run2015B
		// if(isRun2015B)
		// {
		// if(!patUtils::exclusiveDataEventFilter(ev.eventAuxiliary().run(), isMC, isPromptReco ) ) continue;
		// }
		
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
		/*
		bool eTrigger    (
			isMC ? 
			utils::passTriggerPatterns (tr, "HLT_Ele27_eta2p1_WP75_Gsf_v*")
			:
			utils::passTriggerPatterns (tr, "HLT_Ele27_eta2p1_WPLoose_Gsf_v*")
			);
		*/
		bool eTrigger    ( utils::passTriggerPatterns(tr, "HLT_Ele23_WPLoose_Gsf*") );
		bool muTrigger   (
			utils::passTriggerPatterns (tr, "HLT_IsoMu20_v*", "HLT_IsoTkMu20_v*")
			);

		// if(!isMC && muTrigger) mon.fillHisto("nvtx_singlemu_pileup", tags, nGoodPV, 1.);
		// if(!isMC && eTrigger)  mon.fillHisto("nvtx_singlee_pileup",  tags, nGoodPV, 1.);
		
		if(filterOnlySINGLEMU) {                    eTrigger = false; }
		if(filterOnlySINGLEE)  { muTrigger = false;                   }
		
		if (!(eTrigger || muTrigger)) continue;   //ONLY RUN ON THE EVENTS THAT PASS OUR TRIGGERS

		if(debug){
			cout << "Set triggers" << endl;
			}

		// ------------------------------------------------- Apply MET filters
		//if( !isMC && !metFiler.passMetFilter( ev, isPromptReco)) continue;
		

		if(debug){
			cout << "met filters are commented out here" << endl;
			}


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
		// TODO: what is this??
		// there was some wague answer from Pietro.....
		/*
				if(isMC)
					{
						// FIXME: Considering add support for different generators (based on PYTHIA6) for comparison.
						for(size_t igen=0; igen<gen.size(); igen++){
							// FIXME: Should pass to the new status scheme from: https://github.com/cms-sw/cmssw/pull/7791
							//            ////// if(!gen[igen].isHardProcess() && !gen[igen].isPromptFinalState()) continue;
							
							if(gen[igen].status() != 1 &&  gen[igen].status() !=2 && gen[igen].status() !=62 ) continue;
							int absid=abs(gen[igen].pdgId());
							// OK, so taus should be checked as status 2, and quarks as 71 or 23. More testing needed
							//if( absid==15 && hasWasMother(gen[igen]) ) cout << "Event " << iev << ", Particle " << igen << " has " << gen[igen].numberOfDaughters() << " daughters, pdgId " << gen[igen].pdgId() << " and status " << gen[igen].status() << ", mothers " << gen[igen].numberOfMothers() << ", pt " << gen[igen].pt() << ", eta " << gen[igen].eta() << ", phi " << gen[igen].phi() << ". isHardProcess is " << gen[igen].isHardProcess() << ", and isPromptFinalState is " << gen[igen].isPromptFinalState() << endl;
							
							
							//////            if(absid==6 && gen[igen].isHardProcess()){ // particles of the hardest subprocess 22 : intermediate (intended to have preserved mass)
							if(absid==6 && gen[igen].status()==62){ // particles of the hardest subprocess 22 : intermediate (intended to have preserved mass). Josh says 62 (last in chain)
								hasTop=true;

								// FIXME: Top pT reweighting. 13 TeV values not propagated yet, so not using.
								//if(isTTbarMC){
								//  if(gen[igen].get("id") > 0) tPt=gen[igen].pt();
								//  else                        tbarPt=gen[igen].pt();
								//}
							} 
							
							
							//if(!gen[igen].isPromptFinalState() ) continue;
							if( (gen[igen].status() != 1 && gen[igen].status()!= 2 ) || !hasWasMother(gen[igen])) continue;
							
							if((absid==11 || absid==13) && hasLeptonAsDaughter(gen[igen])) cout << "Electron or muon " << igen << " has " << gen[igen].numberOfDaughters() << " daughter which is a lepton." << endl;
							
							if((absid==11 || absid==13) && gen[igen].status()==1)
								{
									ngenLeptonsStatus3++;
									
									if(!hasTauAsMother(gen[igen]))
										ngenLeptonsNonTauSonsStatus3++;
								}
							if(absid==15 && gen[igen].status()==2 )
								{
									ngenTausStatus3++; // This should be summed to ngenLeptonsStatus3 for the dilepton final states, not summed for the single lepton final states.
									//    if(hasLeptonAsDaughter(gen[igen])) cout << "Tau " << igen << " has " << gen[igen].numberOfDaughters() << " daughter which is a lepton." << endl;
								}
							if(absid<=5              ) ngenQuarksStatus3++;
						}
						
						if(debug && (ngenTausStatus3==1 && ngenLeptonsStatus3==1 )  ) cout << "Event: " << iev << ". Leptons: " << ngenLeptonsStatus3 << ". Leptons notaus: " << ngenLeptonsNonTauSonsStatus3 << ". Taus: " << ngenTausStatus3 << ". Quarks: " << ngenQuarksStatus3 << endl;
						
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
						if( 
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
				//      if(tPt>0 && tbarPt>0 && topPtWgt)
				//        {
				//          topPtWgt->computeWeight(tPt,tbarPt);
				//          topPtWgt->getEventWeight(wgtTopPt, wgtTopPtUp, wgtTopPtDown);
				//          wgtTopPtUp /= wgtTopPt;
				//          wgtTopPtDown /= wgtTopPt;
				//        }


					// ------------------------------------ actual particles?

		pat::MuonCollection muons;
		fwlite::Handle<pat::MuonCollection> muonsHandle;
		muonsHandle.getByLabel(ev, "slimmedMuons");
		if(muonsHandle.isValid() ) muons = *muonsHandle;

		pat::ElectronCollection electrons;
		fwlite::Handle<pat::ElectronCollection> electronsHandle;
		electronsHandle.getByLabel(ev, "slimmedElectrons");
		if(electronsHandle.isValid() ) electrons = *electronsHandle;

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
		std::vector<patUtils::GenericLepton> leptons;
		for(size_t l=0; l<electrons.size(); ++l) leptons.push_back(patUtils::GenericLepton (electrons[l] ));
		for(size_t l=0; l<muons.size(); ++l)     leptons.push_back(patUtils::GenericLepton (muons[l]     ));
		std::sort(leptons.begin(), leptons.end(), utils::sort_CandidatesByPt);


				// ---------------------------------- leptons selection
				LorentzVector muDiff(0., 0., 0., 0.);
				std::vector<patUtils::GenericLepton> selLeptons;
				unsigned int nVetoE(0), nVetoMu(0);
				for(size_t ilep=0; ilep<leptons.size (); ++ilep)
					{
						patUtils::GenericLepton& lepton = leptons[ilep];

						bool 
							passKin(true),     passId(true),     passIso(true),
							passVetoKin(true), passVetoId(true), passVetoIso(true);
						
						int lid(lepton.pdgId());
						
						//apply muon corrections
						if(abs(lid) == 13)
						{
							if(muCor)
								{
									TLorentzVector p4(lepton.px(), lepton.py(), lepton.pz(), lepton.energy());
									muCor->applyPtCorrection(p4, lid < 0 ? -1 : 1);
									if(isMC) muCor->applyPtSmearing(p4, lid < 0 ? -1 : 1, false);
									muDiff -= lepton.p4();
									lepton.setP4(LorentzVector(p4.Px(), p4.Py(), p4.Pz(), p4.E()));
									muDiff += lepton.p4();
								}
						}
						
					//no need for charge info any longer
					lid = abs(lid);
					TString lepStr(lid == 13 ? "mu" : "e");
					
					// no need to mess with photon ID // //veto nearby photon (loose electrons are many times photons...)
					// no need to mess with photon ID // double minDRlg(9999.);
					// no need to mess with photon ID // for(size_t ipho=0; ipho<selPhotons.size(); ipho++)
					// no need to mess with photon ID //   minDRlg=TMath::Min(minDRlg,deltaR(leptons[ilep].p4(),selPhotons[ipho].p4()));
					// no need to mess with photon ID // if(minDRlg<0.1) continue;

					// ---------------------------- kinematics
					double leta(fabs(lid==11 ? lepton.el.superCluster()->eta() : lepton.eta()));
					
					// ---------------------- Main lepton kin
					if(lepton.pt() < 30.)                      passKin = false;
					if(leta > 2.1)                                    passKin = false;
					if(lid == 11 && (leta > 1.4442 && leta < 1.5660)) passKin = false; // Crack veto
					
					// ---------------------- Veto lepton kin
					if (lepton.pt () < 20)                      passVetoKin = false;
					if (leta > 2.1)                                    passVetoKin = false;
					if (lid == 11 && (leta > 1.4442 && leta < 1.5660)) passVetoKin = false; // Crack veto
					
					//Cut based identification 
					
					//std::vector<pat::Electron> dummyShit; dummyShit.push_back(leptons[ilep].el);

					// ------------------------- lepton IDs
					passId     = lid == 11 ? patUtils::passId(electronVidMainId, myEvent, lepton.el) : patUtils::passId(lepton.mu, goodPV, patUtils::llvvMuonId::StdTight);
					passVetoId = lid == 11 ? patUtils::passId(electronVidVetoId, myEvent, lepton.el) : patUtils::passId(lepton.mu, goodPV, patUtils::llvvMuonId::StdLoose);

					// ------------------------- lepton isolation
					passIso     = lid == 11 ? true : patUtils::passIso(lepton.mu, patUtils::llvvMuonIso::Tight); // Electron iso is included within the ID
					passVetoIso = lid == 11 ? true : patUtils::passIso(lepton.mu, patUtils::llvvMuonIso::Loose); // Electron iso is included within the ID

					if     (passKin     && passId     && passIso)     selLeptons.push_back(lepton);
					else if(passVetoKin && passVetoId && passVetoIso) lid==11 ? nVetoE++ : nVetoMu++;
					
				}
			std::sort(selLeptons.begin(),   selLeptons.end(),   utils::sort_CandidatesByPt);
			LorentzVector recoMET = met;// FIXME REACTIVATE IT - muDiff;
			
			
			// ------------------------------------------ select the taus
			pat::TauCollection selTaus;
			int ntaus (0);
			for (size_t itau = 0; itau < taus.size(); ++itau)
				{
					pat::Tau& tau = taus[itau];
					if (tau.pt() < 20. || fabs (tau.eta()) > 2.3) continue;
					
					bool overlapWithLepton(false);
					for(int l=0; l<(int)selLeptons.size();++l){
						if(reco::deltaR(tau, selLeptons[l])<0.4){overlapWithLepton=true; break;}
					}
					if(overlapWithLepton) continue;
					
					//      if(!tau.isPFTau()) continue; // Only PFTaus // It should be false for slimmedTaus
					//      if(tau.emFraction() >=2.) continue;
					
					// Discriminators from https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV
					// "The tau passes the discriminator if pat::Tau::tauID("name") returns a value of 0.5 or greater"
					if(tau.tauID("decayModeFindingNewDMs")<0.5) continue; // High pt tau. Otherwise, OldDMs
					// Anyways, the collection of taus from miniAOD should be already afer decayModeFinding cut (the tag - Old or New - is unspecified in the twiki, though).
					// Consequently, there might be a small bias due to events that are cut by the OldDM and would not be cut by the NewDM
					if (tau.tauID ("byMediumCombinedIsolationDeltaBetaCorr3Hits")<0.5) continue; // See whether to us the new byMediumPileupWeightedIsolation3Hits that is available only for dynamic strip reconstruction (default in CMSSW_7_4_14)
					if (tau.tauID ("againstMuonTight3")                          <0.5) continue; // Medium working point not usable. Available values: Loose, Tight
					if (tau.tauID ("againstElectronMediumMVA5")                  <0.5) continue; // Tight working point not usable. Avaiable values: VLoose, Loose, Medium
					
					// Pixel hits cut (will be available out of the box in new MINIAOD production)
					{
						int nChHadPixelHits = 0;
						reco::CandidatePtrVector chCands = tau.signalChargedHadrCands();
						for(reco::CandidatePtrVector::const_iterator iter = chCands.begin(); iter != chCands.end(); iter++){
							pat::PackedCandidate const* packedCand = dynamic_cast<pat::PackedCandidate const*>(iter->get());
							int pixelHits = packedCand->numberOfPixelHits();
							if(pixelHits > nChHadPixelHits) nChHadPixelHits = pixelHits;
						}
						if(nChHadPixelHits==0) continue;
					}
					/////
					
					selTaus.push_back(tau);
					ntaus++;
				}
			std::sort (selTaus.begin(), selTaus.end(), utils::sort_CandidatesByPt);
			
			//
			// ----------------------------------------------- JET/MET ANALYSIS
			//
			if(debug) cout << "Now update Jet Energy Corrections" << endl;
			//add scale/resolution uncertainties and propagate to the MET      
			utils::cmssw::updateJEC(jets,jesCor,totalJESUnc,rho,nGoodPV,isMC);
			if(debug) cout << "Update also MET" << endl;
			std::vector<LorentzVector> newMet=utils::cmssw::getMETvariations(met/*recoMet*/,jets,selLeptons,isMC); // FIXME: Must choose a lepton collection. Perhaps loose leptons?
			met=newMet[utils::cmssw::METvariations::NOMINAL];
			if(debug) cout << "Jet Energy Corrections updated" << endl;
			
			// Select the jets. I need different collections because of tau cleaning, but this is needed only for the single lepton channels, so the tau cleaning is performed later.
			pat::JetCollection
				selJets, selBJets;
			double mindphijmet (9999.);
			for (size_t ijet = 0; ijet < jets.size(); ++ijet)
				{
					pat::Jet& jet = jets[ijet];
					
					if (jet.pt() < 15 || fabs (jet.eta()) > 3.0) continue; // Was 4.7 in eta. Tightened for computing time. 3.0 ensures that we don't cut associations with leptons (0.4 from 2.4)
					
					//mc truth for this jet
					const reco::GenJet * genJet = jet.genJet();
					TString jetType (genJet && genJet->pt() > 0 ? "truejetsid" : "pujetsid");

					//cross-clean with selected leptons and photons
					double minDRlj (9999.), minDRlg (9999.), minDRljSingleLep(9999.);

					for (size_t ilep = 0; ilep < selLeptons.size(); ilep++)
						minDRlj = TMath::Min(minDRlj, reco::deltaR (jet, selLeptons[ilep]));
					// don't want to mess with photon ID // for(size_t ipho=0; ipho<selPhotons.size(); ipho++)
					// don't want to mess with photon ID //   minDRlg = TMath::Min( minDRlg, deltaR(jets[ijet],selPhotons[ipho]) );
					//          if (minDRlj < 0.4 /*|| minDRlg<0.4 */ ) continue;
					
					for (size_t ilep = 0; ilep < selLeptons.size(); ilep++)
						minDRljSingleLep = TMath::Min(minDRljSingleLep, reco::deltaR (jet, selLeptons[ilep]));
					
					//jet id
					bool passPFloose = passPFJetID("Loose", jet); 
					// FIXME: check when pileup ID will come out
					//if (jets[ijet].pt() > 30)
					//  {
					//    mon.fillHisto (jetType, "", fabs (jets[ijet].eta()), 0);
					//    if (passPFloose)                        mon.fillHisto (jetType, "", fabs (jets[ijet].eta()), 1);
					//    if (passLooseSimplePuId)                mon.fillHisto (jetType, "", fabs (jets[ijet].eta()), 2);
					//    if (passPFloose && passLooseSimplePuId) mon.fillHisto (jetType, "", fabs (jets[ijet].eta()), 3);
					//  }
					if (!passPFloose || jet.pt() <30. || fabs(jet.eta()) > 2.5) continue;
					if (minDRlj < 0.4) continue;
					
					selJets.push_back(jet);
					
					double dphijmet = fabs (deltaPhi (met.phi(), jet.phi()));
					if (dphijmet < mindphijmet) mindphijmet = dphijmet;
					bool hasCSVtag(jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > btagMedium);
					
					if (isMC)
						{
							int flavId = jet.partonFlavour();
							if      (abs(flavId)==5) btsfutil.modifyBTagsWithSF(hasCSVtag, sfb,   beff);
							else if (abs(flavId)==4) btsfutil.modifyBTagsWithSF(hasCSVtag, sfb/5, beff);
							else                     btsfutil.modifyBTagsWithSF(hasCSVtag, sfl,   leff);
						}
					
					if(!hasCSVtag) continue;
					selBJets.push_back(jet);
				}

			std::sort (selJets.begin(),  selJets.end(),  utils::sort_CandidatesByPt);
			std::sort (selBJets.begin(), selBJets.end(), utils::sort_CandidatesByPt);
			
			//
			// -------------------------------------------------- ASSIGN CHANNEL
			//
			std::vector < TString > chTags; chTags.clear();
			int 
				dilId (1),
				slepId(0);
			LorentzVector dileptonSystem (0, 0, 0, 0);
			if(selLeptons.size()>=2)
				{
					for (size_t ilep = 0; ilep < 2; ilep++)
						{
							dilId *= selLeptons[ilep].pdgId();
							int id(abs (selLeptons[ilep].pdgId()));
							dileptonSystem += selLeptons[ilep].p4();
						}
				}
			
			if(selLeptons.size()>0)
				slepId=selLeptons[0].pdgId();
			
			// Event classification. Single lepton triggers are used for offline selection of dilepton events. The "else if"s guarantee orthogonality
			bool 
				isSingleMu(false),
				isSingleE(false),
				isDoubleMu(false),
				isDoubleE(false),
				isEMu(false);
			int multiChannel(0);
			if      (abs(slepId)==13 && muTrigger && nVetoE==0 && nVetoMu==0){ isSingleMu = true; multiChannel++; chTags.push_back("singlemu");}
			else if (abs(slepId)==11 && eTrigger  && nVetoE==0 && nVetoMu==0){ isSingleE  = true; multiChannel++; chTags.push_back("singlee");}
			else if (abs(dilId)==121 && eTrigger                            ){ isDoubleE  = true; multiChannel++; chTags.push_back("ee");}
			else if (abs(dilId)==169 && muTrigger                           ){ isDoubleMu = true; multiChannel++; chTags.push_back("mumu");}
			else if (abs(dilId)==143 && muTrigger                           ){ isEMu      = true; multiChannel++; chTags.push_back("emu");}
			else if (abs(dilId)==143 && eTrigger  && !muTrigger             ){ isEMu      = true; multiChannel++; chTags.push_back("emu");} // Pick up the largest number of emu events possible, maintaining orthogonality
			
			// keep in mind the eventCategory thingy for more refined categorization // TString evCat=eventCategoryInst.GetCategory(selJets,dileptonSystem);
			//std::vector < TString > tags (1, "all");
			for (size_t ich = 0; ich < chTags.size(); ich++)
				{
					tags.push_back (chTags[ich]);
					//tags.push_back( chTags[ich]+evCat );
				}
			if(multiChannel>1) nMultiChannel++;

		if(debug){
			cout << "channel is defined, running the event selection" << endl;
		}

			// Dilepton full analysis
			// No dilepton analysis
			//if( isDoubleE || isEMu || isDoubleMu){ continue; }

			// ------------------------------------------ Single lepton full analysis
			//if(tags[1] == "singlemu" || tags[1] == "singlee")
		if(isSingleMu || isSingleE){
			singlelep_ttbar_preselectedevents->Fill(1);

			// ---------------------------- Clean jet collection from selected taus
			pat::JetCollection
			selSingleLepJets, selSingleLepBJets;
			for (size_t ijet = 0; ijet < selJets.size(); ++ijet)
				{
				pat::Jet jet = selJets[ijet];

				double minDRtj(9999.);
				for(size_t itau=0; itau<selTaus.size(); ++itau)
					{
					minDRtj = TMath::Min(minDRtj, reco::deltaR(jet, selTaus[itau]));
					}
				if(minDRtj>0.4) selSingleLepJets.push_back(jet);

				bool hasCSVtag = (jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > btagMedium);
				if (isMC)
					{
					int flavId = jets[ijet].partonFlavour();
					if      (abs (flavId) == 5) btsfutil.modifyBTagsWithSF(hasCSVtag, sfb,   beff);
					else if (abs (flavId) == 4) btsfutil.modifyBTagsWithSF(hasCSVtag, sfb/5, beff);
					else                        btsfutil.modifyBTagsWithSF(hasCSVtag, sfl,   leff);
					}

				if(!hasCSVtag) continue;
				if(minDRtj>0.4) selSingleLepBJets.push_back(jets[ijet]);
				}

			std::sort(selSingleLepJets.begin(),  selSingleLepJets.end(),  utils::sort_CandidatesByPt);
			std::sort(selSingleLepBJets.begin(), selSingleLepBJets.end(), utils::sort_CandidatesByPt);


			// mon.fillHisto("nvtx_pileup", tags, nGoodPV, weight);

			if(selLeptons.size()!=1 || nGoodPV==0) continue; // Veto requirement alredy applied during the event categoriziation
			//int id (abs (selLeptons[0].pdgId()));
			//weight *= isMC ? lepEff.getLeptonEfficiency(selLeptons[0].pt(), selLeptons[0].eta(), id, id == 11 ? "loose" : "tight").first : 1.0;        

			// Event selection booleans
			bool passJetSelection(selSingleLepJets.size()>1); // 2 jets
			bool passMetSelection(met.pt()>40.); // MET > 40
			bool passBtagsSelection(selSingleLepBJets.size()>0); // 1 b jet
			bool passTauSelection(selTaus.size()==1); // only 1 tau
			bool passOS(selTaus.size()>0 ? selLeptons[0].pdgId() * selTaus[0].pdgId() < 0 : 0); // Oposite sign

			if (passJetSelection)
				{
				if(isSingleMu) singlelep_ttbar_selected2_mu_events->Fill(1);
				else if (isSingleE) singlelep_ttbar_selected2_el_events->Fill(1);
				}


			// Mara's selection booleans
			//bool passMaraJetSelection(selSingleLepJets.size()>3); // 4 jets FIXME: why 4 jets?
			//bool passMaraBtagsSelection(selSingleLepBJets.size()>1); // 2 b-tag
			//bool passMaraLeptonSelection( selLeptons.size()>0 ); // 1 lepton -- should be included by default at this stage

			// common-selection
			if(passJetSelection && passBtagsSelection) // 2 jets, 1 b jet, 1 isolated lepton
				{
				if(isSingleMu) singlelep_ttbar_selected_mu_events->Fill(1);
				else if (isSingleE) singlelep_ttbar_selected_el_events->Fill(1);
				crossel_sum_weights_raw += rawWeight;
				crossel_sum_weights += weight;
				/*
				fprintf(csv_out, "crossel:%d,%d,%g,%g,%d,", num_inters, nGoodPV, rawWeight, weight, isSingleE);
				pb.SetPxPyPzE( selSingleLepBJets[0].px(), selSingleLepBJets[0].py(), selSingleLepBJets[0].pz(), selSingleLepBJets[0].pt()); // 
				pbb.SetPxPyPzE( selSingleLepJets[1].px(), selSingleLepJets[1].py(), selSingleLepJets[1].pz(), selSingleLepJets[1].pt()); // or take another B???
				pl.SetPxPyPzE( selLeptons[0].px(), selLeptons[0].py(), selLeptons[0].pz(), selLeptons[0].pt());
				plb.SetPxPyPzE( selLeptons[1].px(), selLeptons[1].py(), selLeptons[1].pz(), selLeptons[1].pt());
				// no kino1 for now
				//prest = pb + pbb + pl + plb;
				//fprintf(csv_out, "kino1:\n");
				//fprintf(csv_out, "%g, %g, %g, %g,", pl.E(), plb.E(), pb.E(), pbb.E());
				//fprintf(csv_out, "%g, %g, %g,", (prest.X()*prest.X() + prest.Y()*prest.Y() + prest.Z()*prest.Z()),
						//(prest.X()*prest.X() + prest.Y()*prest.Y()), met.pt());
				//fprintf(csv_out, "%g, %g, %g,\n", prest*(pl+pb), pl*pb, plb*pbb);
				//fprintf(csv_out, "%g, %g, %g\n", (pl+pb).Angle(prest.Vect()), pl.Angle(pb.Vect()), plb.Angle(pbb.Vect()));
				//fprintf(csv_out, "kino2:\n");
				fprintf(csv_out, "%g,%g,%g,%g,", pl.X(), pl.Y(), pl.Z(), pl.E());
				//fprintf(csv_out, "%g,%g,%g,%g,", plb.X(), plb.Y(), plb.Z(), plb.E());
				fprintf(csv_out, "%g,%g,%g,%g,", pb.X(), pb.Y(), pb.Z(), pb.E());
				//fprintf(csv_out, "%g,%g,%g,%g\n", pbb.X(), pbb.Y(), pbb.Z(), pbb.E());
				fprintf(csv_out, "%g,%g,%g,%g,", selSingleLepJets[0].px(), selSingleLepJets[0].py(), selSingleLepJets[0].pz(), selSingleLepJets[0].pt() );
				fprintf(csv_out, "%g,%g,%g,%g\n", selSingleLepJets[1].px(), selSingleLepJets[1].py(), selSingleLepJets[1].pz(), selSingleLepJets[1].pt() );
				*/
				}


			// oursel: 2 jets, 1 b, 1 iso lepton, 1 tau
			if(passJetSelection && passMetSelection && passBtagsSelection && passTauSelection && passOS )
				{
				if(isSingleMu) singlelep_ttbar_selected_mu_events->Fill(1);
				else if (isSingleE) singlelep_ttbar_selected_el_events->Fill(1);
				fprintf(csv_out, "oursel:%d,%d,%g,%g,%d,", num_inters, nGoodPV, rawWeight, weight, isSingleE);
				oursel_sum_weights_raw += rawWeight;
				oursel_sum_weights += weight;
				pb.SetPxPyPzE( selSingleLepBJets[0].px(), selSingleLepBJets[0].py(), selSingleLepBJets[0].pz(), selSingleLepBJets[0].pt()); // 
				//pbb.SetPxPyPzE( selSingleLepJets[1].px(), selSingleLepJets[1].py(), selSingleLepJets[1].pz(), selSingleLepJets[1].pt()); // or take another B???
				pl.SetPxPyPzE( selLeptons[0].px(), selLeptons[0].py(), selLeptons[0].pz(), selLeptons[0].pt());
				//plb.SetPxPyPzE( selLeptons[1].px(), selLeptons[1].py(), selLeptons[1].pz(), selLeptons[1].pt());
				//prest = pb + pbb + pl + plb;
				//fprintf(csv_out, "kino1:\n");
				//fprintf(csv_out, "%g, %g, %g, %g,\n", pl.E(), plb.E(), pb.E(), pbb.E());
				//fprintf(csv_out, "%g, %g, %g,\n", (prest.X()*prest.X() + prest.Y()*prest.Y() + prest.Z()*prest.Z()),
						//(prest.X()*prest.X() + prest.Y()*prest.Y()), met.pt());
				//fprintf(csv_out, "%g, %g, %g,\n", prest*(pl+pb), pl*pb, plb*pbb);
				//fprintf(csv_out, "%g, %g, %g\n", (pl+pb).Angle(prest.Vect()), pl.Angle(pb.Vect()), plb.Angle(pbb.Vect()));
				//fprintf(csv_out, "kino2:\n");
				fprintf(csv_out, "%g,%g,%g,%g,",  pl.X(), pl.Y(), pl.Z(), pl.E());
				fprintf(csv_out, "%g,%g,%g,%g,", selTaus[0].px(), selTaus[0].py(), selTaus[0].pz(), selTaus[0].pt() );
				//fprintf(csv_out, "%g,%g,%g,%g,", plb.X(), plb.Y(), plb.Z(), plb.E());
				fprintf(csv_out, "%g,%g,%g,%g,",  pb.X(), pb.Y(), pb.Z(), pb.E());
				fprintf(csv_out, "%g,%g,%g,%g,", selSingleLepJets[0].px(), selSingleLepJets[0].py(), selSingleLepJets[0].pz(), selSingleLepJets[0].pt() );
				fprintf(csv_out, "%g,%g,%g,%g\n", selSingleLepJets[1].px(), selSingleLepJets[1].py(), selSingleLepJets[1].pz(), selSingleLepJets[1].pt() );
				//fprintf(csv_out, "%g,%g,%g,%g\n", pbb.X(), pbb.Y(), pbb.Z(), pbb.E());
				}

			/*
			pb.SetVect( selSingleLepBJets[0].momentum() ); // 
			pb.E = selSingleLepBJets[0].correctedP4();
			pbb.SetVect( selSingleLepJets[1].momentum() ); // or take another B???
			pbb.E = selSingleLepJets[1].correctedP4();
			pl = selLeptons[0];
			plb = selLeptons[1];
			prest = pb + pbb + pl + plb;
			//kino->Fill( el, elb, eb, ebb, prest, prestt, met, thetascal, phiscal, phibscal, theta, phi, phib );
			// pl - pb, plb - pbb
			kino->Fill( pl.E(), plb.E(), pb.E(), pbb.E(), (prest.X()*prest.X() + prest.Y()*prest.Y() + prest.Z()*prest.Z()),
								 (prest.X()*prest.X() + prest.Y()*prest.Y()), met.pt(),
									prest*(pl+pb), pl*pb, plb*pbb,
								 (pl+pb).Angle(prest.Vect()), pl.Angle(pb.Vect()), plb.Angle(pbb.Vect()) );

					 // reverse: pl - pbb
					 kino->Fill( plb.E(), pl.E(), pb.E(), pbb.E(), (prest.X()*prest.X() + prest.Y()*prest.Y() + prest.Z()*prest.Z()),
											(prest.X()*prest.X() + prest.Y()*prest.Y()), met.pt(),
											 prest*(plb+pb), plb*pb, pl*pbb,
											(plb+pb).Angle(prest.Vect()), plb.Angle(pb.Vect()), pl.Angle(pbb.Vect()) );
					 //kino->Fill( pl, el, plb, elb, pb, eb, pbb, ebb, prest, prestt, evmet, theta, phi, thetab, phib );
					 //// and the other combination
					 //kino->Fill( pl, el, plb, elb, pb, eb, pbb, ebb, prest, prestt, evmet, theta, phi, thetab, phib );
			*/

			// Mara's selection booleans
			bool passMaraJetSelection(selSingleLepJets.size()>3); // 4 jets
			bool passMaraBtagsSelection(selSingleLepBJets.size()>1); // 2 b-tag
			bool passMaraLeptonSelection( selLeptons.size()>0 ); // 1 lepton

			if(passMaraJetSelection && passMaraBtagsSelection && passMaraLeptonSelection)
				{
				if(isSingleMu) singlelep_ttbar_maraselected_mu_events->Fill(1);
				else if (isSingleE) singlelep_ttbar_maraselected_el_events->Fill(1);
				marasel_sum_weights_raw += rawWeight;
				marasel_sum_weights += weight;
				fprintf(csv_out, "marasel:%d,%d,%g,%g,%d,", num_inters, nGoodPV, rawWeight, weight, isSingleE);
				// TODO: print out separately b-jets and all other jets?
				pb.SetPxPyPzE(  selSingleLepBJets[0].px(), selSingleLepBJets[0].py(), selSingleLepBJets[0].pz(), selSingleLepBJets[0].pt()); // 
				pbb.SetPxPyPzE( selSingleLepBJets[1].px(), selSingleLepBJets[1].py(), selSingleLepBJets[1].pz(), selSingleLepBJets[1].pt());
				pl.SetPxPyPzE(  selLeptons[0].px(), selLeptons[0].py(), selLeptons[0].pz(), selLeptons[0].pt());
				//plb.SetPxPyPzE( selLeptons[1].px(), selLeptons[1].py(), selLeptons[1].pz(), selLeptons[1].pt());
				//fprintf(csv_out, "kino2:\n");
				fprintf(csv_out, "%g,%g,%g,%g,", pl.X(), pl.Y(), pl.Z(), pl.E());
				//fprintf(csv_out, "%g,%g,%g,%g,", plb.X(), plb.Y(), plb.Z(), plb.E());
				fprintf(csv_out, "%g,%g,%g,%g,", pb.X(), pb.Y(), pb.Z(), pb.E());
				fprintf(csv_out, "%g,%g,%g,%g,", pbb.X(), pbb.Y(), pbb.Z(), pbb.E());
				fprintf(csv_out, "%g,%g,%g,%g,", selSingleLepJets[0].px(), selSingleLepJets[0].py(), selSingleLepJets[0].pz(), selSingleLepJets[0].pt() );
				fprintf(csv_out, "%g,%g,%g,%g,", selSingleLepJets[1].px(), selSingleLepJets[1].py(), selSingleLepJets[1].pz(), selSingleLepJets[1].pt() );
				fprintf(csv_out, "%g,%g,%g,%g,", selSingleLepJets[2].px(), selSingleLepJets[2].py(), selSingleLepJets[2].pz(), selSingleLepJets[2].pt() );
				fprintf(csv_out, "%g,%g,%g,%g\n", selSingleLepJets[3].px(), selSingleLepJets[3].py(), selSingleLepJets[3].pz(), selSingleLepJets[3].pt() );
				}

				/* // old crap with smartmon:
				// Setting up control categories and fill up event flow histo
				std::vector < TString > ctrlCats;
				ctrlCats.clear ();
																																																			{ ctrlCats.push_back ("step1"); mon.fillHisto("xseceventflowslep", tags, 0, weight); mon.fillHisto("chhiggseventflowslep", tags, 0, weight); }
				if(passJetSelection   )                                                                       { ctrlCats.push_back ("step2"); mon.fillHisto("xseceventflowslep", tags, 1, weight); mon.fillHisto("chhiggseventflowslep", tags, 1, weight);
					 if(isSingleMu) singlelep_ttbar_selected2_mu_events->Fill(1);
					 else if (isSingleE) singlelep_ttbar_selected2_el_events->Fill(1);
				}
				if(passJetSelection && passMetSelection )                                                     { ctrlCats.push_back ("step3"); mon.fillHisto("xseceventflowslep", tags, 2, weight); mon.fillHisto("chhiggseventflowslep", tags, 2, weight); }
				if(passJetSelection && passMetSelection && passBtagsSelection )                               { ctrlCats.push_back ("step4"); mon.fillHisto("xseceventflowslep", tags, 3, weight); mon.fillHisto("chhiggseventflowslep", tags, 3, weight); }
				if(passJetSelection && passMetSelection && passBtagsSelection && passTauSelection )           { ctrlCats.push_back ("step5"); mon.fillHisto("xseceventflowslep", tags, 4, weight); mon.fillHisto("chhiggseventflowslep", tags, 4, weight); }
				if(passJetSelection && passMetSelection && passBtagsSelection && passTauSelection && passOS ) { ctrlCats.push_back ("step6"); mon.fillHisto("xseceventflowslep", tags, 5, weight); mon.fillHisto("chhiggseventflowslep", tags, 5, weight);
					 if(isSingleMu) singlelep_ttbar_selected_mu_events->Fill(1);
					 else if (isSingleE) singlelep_ttbar_selected_el_events->Fill(1);
				}
				

				bool passBtagsSelection_0(selSingleLepBJets.size()==0);
				bool passBtagsSelection_1(selSingleLepBJets.size()==1);
				bool passBtagsSelection_2(selSingleLepBJets.size()>1);

																																										{ ctrlCats.push_back("altstep1"); mon.fillHisto("xsecalteventflowslep", tags, 0, weight); mon.fillHisto("chhiggsalteventflowslep", tags, 0, weight); }
				if(passMetSelection)                                                        { ctrlCats.push_back("altstep2"); mon.fillHisto("xsecalteventflowslep", tags, 1, weight); mon.fillHisto("chhiggsalteventflowslep", tags, 1, weight); }
				if(passMetSelection && passTauSelection)                                    { ctrlCats.push_back("altstep3"); mon.fillHisto("xsecalteventflowslep", tags, 2, weight); mon.fillHisto("chhiggsalteventflowslep", tags, 2, weight); }
				if(passMetSelection && passTauSelection && passOS)                          { ctrlCats.push_back("altstep4"); mon.fillHisto("xsecalteventflowslep", tags, 3, weight); mon.fillHisto("chhiggsalteventflowslep", tags, 3, weight); }
				if(passMetSelection && passTauSelection && passOS && passBtagsSelection_0)  { ctrlCats.push_back("altstep5"); mon.fillHisto("xsecalteventflowslep", tags, 4, weight); mon.fillHisto("chhiggsalteventflowslep", tags, 4, weight); }
				if(passMetSelection && passTauSelection && passOS && passBtagsSelection_1)  { ctrlCats.push_back("altstep6"); mon.fillHisto("xsecalteventflowslep", tags, 5, weight); mon.fillHisto("chhiggsalteventflowslep", tags, 5, weight); }
				if(passMetSelection && passTauSelection && passOS && passBtagsSelection_2)  { ctrlCats.push_back("altstep7"); mon.fillHisto("xsecalteventflowslep", tags, 6, weight); mon.fillHisto("chhiggsalteventflowslep", tags, 6, weight); }

				// Fill the control plots
				for(size_t k=0; k<ctrlCats.size(); ++k){
					
					TString icat(ctrlCats[k]);
					mon.fillHisto(icat+"nvtxraw",    tags, nGoodPV,                            rawWeight);
					mon.fillHisto(icat+"nvtx",       tags, nGoodPV,                            weight   );
					mon.fillHisto(icat+"rho",        tags, rho,                                weight   );
					mon.fillHisto(icat+"leadpt",     tags, selLeptons[0].pt(),                 weight   );
					mon.fillHisto(icat+"leadleptonpt",     tags, selLeptons[0].pt(),                 weight   );
					mon.fillHisto(icat+"trailerpt",  tags, selLeptons[1].pt(),                 weight   );
					mon.fillHisto(icat+"trailerleptonpt",  tags, selLeptons[1].pt(),                 weight   );
					mon.fillHisto(icat+"leadleptoneta",    tags, fabs(selLeptons[0].eta()),          weight   );
					mon.fillHisto(icat+"trailerleptoneta", tags, fabs(selLeptons[1].eta()),          weight   );
					mon.fillHisto(icat+"ntaus",      tags, ntaus,                              weight   );
					mon.fillHisto(icat+"met",        tags, met.pt(),                           weight   );
					mon.fillHisto(icat+"recomet",    tags, recoMET.pt(),                       weight   );
					if(selSingleLepJets.size()>0){
						mon.fillHisto(icat+"leadjetpt",      tags, selSingleLepJets[0].pt(),         weight);
						//mon.fillHisto(icat+"trailerpt",   tags, selLeptons[1].pt(),         weight);
						mon.fillHisto(icat+"leadjeteta",     tags, fabs (selSingleLepJets[0].eta()), weight);
						//mon.fillHisto(icat+"trailereta",  tags, fabs (selLeptons[1].eta()), weight);
					}
					if(ntaus > 0){
						mon.fillHisto (icat+"tauleadpt", tags, selTaus[0].pt(),             weight);
						mon.fillHisto (icat+"tauleadeta", tags, selTaus[0].eta(),             weight);
					}
				 
					
					mon.fillHisto(icat+"nbtags", tags, selSingleLepBJets.size(), weight);
					mon.fillHisto(icat+"njets",  tags, selSingleLepJets.size(), weight);
					// dilepton only           mon.fillHisto (icat+"zmass", tags, dileptonSystem.mass(),           weight);
					// dilepton only           mon.fillHisto (icat+"zy",    tags, fabs(dileptonSystem.Rapidity()), weight);
					// dilepton only           mon.fillHisto (icat+"zpt",   tags, dileptonSystem.pt(),             weight);
					// dilepton only           //these two are used to reweight photon -> Z, the 3rd is a control
					// dilepton only           mon.fillHisto (icat+"qt",    tags, dileptonSystem.pt(),             weight, true);
					// dilepton only           ///     mon.fillHisto("qtraw",    tags, dileptonSystem.pt(),weight/triggerPrescale,true); 
					
					for (size_t ijet = 0; ijet < selSingleLepJets.size(); ijet++)
						{
							if (selSingleLepJets[ijet].pt() < 30 || fabs (selSingleLepJets[ijet].eta()) > 2.5) continue;
							
							double csv (selSingleLepJets[ijet].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
							mon.fillHisto ("csv", tags, csv, weight);
							if (!isMC) continue;
							int flavId = selSingleLepJets[ijet].partonFlavour();
							TString jetFlav ("others");
							if (abs (flavId) == 5) jetFlav = "b";
							else if (abs (flavId) == 4) jetFlav = "c";
							mon.fillHisto ("csv" + jetFlav, tags, csv, weight);
						}
				}
				//
				// HISTOS FOR STATISTICAL ANALYSIS (include systematic variations)
				//
				//Fill histogram for posterior optimization, or for control regions
				
				if(passTauSelection && passOS)
					{
						for (size_t ivar = 0; ivar < nSystVars; ivar++)
							{
								TString var(systVars[ivar]);
								
								double iweight = weight;       //nominal
								
								//energy scale/resolution
								bool varyJesUp    (systVars[ivar] == "_jesup"   );
								bool varyJesDown  (systVars[ivar] == "_jesdown" );
								bool varyJerUp    (systVars[ivar] == "_jerup"   );
								bool varyJerDown  (systVars[ivar] == "_jerdown" );
								bool varyUmetUp   (systVars[ivar] == "_umetup"  );
								bool varyUmetDown (systVars[ivar] == "_umetdown");
								bool varyLesUp    (systVars[ivar] == "_lesup"   );
								bool varyLesDown  (systVars[ivar] == "_lesdown" );
								
								//pileup variations
								if (systVars[ivar] == "_puup")   iweight *= TotalWeight_plus;
								if (systVars[ivar] == "_pudown") iweight *= TotalWeight_minus;

								//btag
								bool varyBtagUp (systVars[ivar] == "_btagup");
								bool varyBtagDown (systVars[ivar] == "_btagdown");
								
								//Here were the Q^2 variations on VV pT spectum
								
								//recompute MET/MT if JES/JER was varied
								LorentzVector newMET = mets[0].p4();
								
								if(varyJesUp)    newMET = mets[0].shiftedP4(pat::MET::METUncertainty::JetEnUp);
								if(varyJesDown)  newMET = mets[0].shiftedP4(pat::MET::METUncertainty::JetEnDown);
								if(varyJerUp)    newMET = mets[0].shiftedP4(pat::MET::METUncertainty::JetResUp);
								if(varyJerDown)  newMET = mets[0].shiftedP4(pat::MET::METUncertainty::JetResDown);
								if(varyUmetUp)   newMET = mets[0].shiftedP4(pat::MET::METUncertainty::UnclusteredEnUp);
								if(varyUmetDown) newMET = mets[0].shiftedP4(pat::MET::METUncertainty::UnclusteredEnDown);
								//if(varyLesUp)    newMET = met[utils::cmssw::LESUP]; //FIXME  must vary all leptons separately: MuonEnUp/MuonEnDown/ElectronEnUp/ElectronEnDown/TauEnUp/TauEnDown
								//if(varyLesDown)  newMET = met[utils::cmssw::LESDOWN];
								
								pat::JetCollection finalSelSingleLepJets;
								pat::JetCollection finalSelSingleLepBJets;
								bool passLocalBveto (true);///passBtags);
								for (size_t ijet = 0; ijet < jets.size(); ijet++)
									{
										pat::Jet jet = jets[ijet];
										double eta = jet.eta();
										double pt = jet.pt();
										if(isMC)
											{
												std::vector<double> varPt = utils::cmssw::smearJES(pt, eta, totalJESUnc);
												if(varyJesUp)   pt = varPt[0];
												if(varyJesDown) pt = varPt[1];
												//  smearJER(float pt, float eta, float genPt)
												//  float newJERSF(1.0);
												//if(isMC)
												//  {
												//    const data::PhysicsObject_t &genJet=jets[ijet].getObject("genJet");
												//    std::vector<float> smearJER=utils::cmssw::smearJER(jets[ijet].pt(),jets[ijet].eta(),genJet.pt());
												//    newJERSF=smearJER[0]/jets[ijet].pt();
												//    rawJet *= newJERSF;
												// if(varyJerUp)    pt=jets[ijet].getVal("jerup");
												// if(varyJerDown)  pt=jets[ijet].getVal("jerdown");
											}
										
										if (pt < 30 || fabs(eta) > 2.5) continue;
										bool passPFloose = passPFJetID("Loose", jet); 
										if (!passPFloose) continue;
										
										//cross-clean with selected leptons and photons
										double minDRlj (9999.), minDRlg (9999.);
										for (size_t ilep = 0; ilep < selLeptons.size(); ilep++)
											minDRlj = TMath::Min (minDRlj, reco::deltaR (jet.p4(), selLeptons[ilep].p4()));
										// don't want to mess with photon ID // for(size_t ipho=0; ipho<selPhotons.size(); ipho++)
										// don't want to mess with photon ID //   minDRlg = TMath::Min( minDRlg, deltaR(jets[ijet].p4(),selPhotons[ipho].p4()) );
										double minDRtj(9999.);
										for(size_t itau=0; itau<selTaus.size(); ++itau)
											{
												minDRtj = TMath::Min(minDRtj, reco::deltaR(jet, selTaus[itau]));
											}
										if (minDRlj < 0.4 || minDRtj<0.4 ) continue;
										
										finalSelSingleLepJets.push_back(jet);
										
										int flavId(jet.partonFlavour());
										bool hasCSVtag(jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > btagMedium);

				if (varyBtagUp)
					{
					if (abs (flavId) == 5)      btsfutil.modifyBTagsWithSF(hasCSVtag, sfb + sfbunc,     beff);
					else if (abs (flavId) == 4) btsfutil.modifyBTagsWithSF(hasCSVtag, sfb/5 + 2*sfbunc, beff);
					else                        btsfutil.modifyBTagsWithSF(hasCSVtag, sfl + sflunc,     leff);
					}
				else if (varyBtagDown)
					{
					if (abs (flavId) == 5)      btsfutil.modifyBTagsWithSF(hasCSVtag, sfb - sfbunc,     beff);
					else if (abs (flavId) == 4) btsfutil.modifyBTagsWithSF(hasCSVtag, sfb/5 - 2*sfbunc, beff);
					else                        btsfutil.modifyBTagsWithSF(hasCSVtag, sfl - sflunc,     leff);
					}
				if(hasCSVtag) finalSelSingleLepBJets.push_back(jet);
				}
				std::sort(finalSelSingleLepJets.begin(),  finalSelSingleLepJets.end(),  utils::sort_CandidatesByPt);
				std::sort(finalSelSingleLepBJets.begin(), finalSelSingleLepBJets.end(), utils::sort_CandidatesByPt);

				bool passFinalJetSelection(finalSelSingleLepJets.size()>1);
				bool passFinalMetSelection(newMET.pt()>40.);
				bool passFinalBtagsSelection(finalSelSingleLepBJets.size()>0);

				if(!passFinalJetSelection || !passFinalMetSelection || !passFinalBtagsSelection) continue;
				// Here fill stat plots
				pat::Tau & tau = selTaus[0];
				reco::CandidatePtr leadChargedHadron = tau.leadChargedHadrCand();
				double tauR(leadChargedHadron->p() / tau.energy());  // Sic. It is momentum, not transverse momentum
				double tauY(2*leadChargedHadron->pt()/tau.et() - 1);
				//                Y= [ p_T^{trk} - (E_T - p_T^{trk} ]/E_T  = 2p_T^{trk}/E_T  - 1 . Which is practically Y = 2R' -1


				LorentzVector mutauSystem (0, 0, 0, 0);
				mutauSystem += selLeptons[0].p4();
				mutauSystem += tau.p4();
								
				mon.fillHisto("finalnbjets"         +var, tags, finalSelSingleLepBJets.size(), iweight);
				mon.fillHisto("finaltaur"           +var, tags, tauR, iweight);
				mon.fillHisto("finaltaupolarization"+var, tags, tauY, iweight);
				mon.fillHisto("finaldphilepmet"     +var, tags, fabs(deltaPhi(newMET.phi(), selLeptons[0].phi())), iweight);
				mon.fillHisto("finaldphitaumet"     +var, tags, fabs(deltaPhi(newMET.phi(), selTaus[0].phi())), iweight);
				mon.fillHisto("finaldphileptau"     +var, tags, fabs(deltaPhi(selLeptons[0].phi(), selTaus[0].phi())), iweight);
				mon.fillHisto("finaltaupt"          +var, tags, selTaus[0].pt(), iweight);
				mon.fillHisto("finalmutaumass"      +var, tags, mutauSystem.mass(), iweight);

				if(saveSummaryTree)
					{
					TDirectory* cwd = gDirectory;
					summaryFile->cd();
					summaryTree->Fill();
					cwd->cd();
					}
				}
			} // End stat analysis

			*/
			} // End single lepton full analysis

		if(debug){
			cout << "Finished processing first event in the first file, exiting" << endl;
			//return 0;
			break;
			}

		} // End single file event loop
	//fprintf(csv_out, "In the file, num_events, sum_weights_raw, sum_weights: %d, %g, %g\n", iev, sum_weights_raw, sum_weights);
	fprintf(csv_out, "acceptances:");
	fprintf(csv_out, "%s,", urls[f].c_str());
	fprintf(csv_out, "%d,%g,%g,", iev, sum_weights_raw, sum_weights);
	fprintf(csv_out, "%g,%g,",  crossel_sum_weights_raw, crossel_sum_weights);
	fprintf(csv_out, "%g,%g,",  oursel_sum_weights_raw, oursel_sum_weights);
	fprintf(csv_out, "%g,%g\n", marasel_sum_weights_raw, marasel_sum_weights);

	fprintf(csv_out, "negative_events_nvtx:");
        for (int i=0; i<100; i++)
                { fprintf(csv_out, "%d,", negative_event_nvtx[i]); }
	fprintf(csv_out, "\n");
	fprintf(csv_out, "positive_events_nvtx:");
        for (int i=0; i<100; i++)
                { fprintf(csv_out, "%d,", positive_event_nvtx[i]); }
	fprintf(csv_out, "\n");
	printf("\n");

	delete file;
	} // End loop on files

fclose(csv_out);

if(saveSummaryTree)
	{
	TDirectory* cwd = gDirectory;
	summaryFile->cd();
	summaryTree->Write();
	summaryFile->Close();
	delete summaryFile;
	cwd->cd();
	}


if(nMultiChannel>0) cout << "Warning! There were " << nMultiChannel << " multi-channel events out of " << totalEntries << " events!" << endl;
printf ("\n");

//##############################################
//########     SAVING HISTO TO FILE     ########
//##############################################
//save control plots to file
printf ("Results save in %s\n", outUrl.Data());

//save all to the file
TFile *ofile = TFile::Open (outUrl + ".root", "recreate");
// mon.Write();

singlelep_ttbar_initialevents->Write();
singlelep_ttbar_preselectedevents->Write();
singlelep_ttbar_selected_mu_events->Write();
singlelep_ttbar_selected_el_events->Write();
singlelep_ttbar_selected2_mu_events->Write();
singlelep_ttbar_selected2_el_events->Write();

singlelep_ttbar_maraselected_mu_events->Write();
singlelep_ttbar_maraselected_el_events->Write();

ofile->Close();

if (outTxtFile) fclose (outTxtFile);

// Now that everything is done, dump the list of lumiBlock that we processed in this job
if(!isMC){
	// FIXME: when lumi certificate is ready for rereco data, check that these work
	goodLumiFilter.FindLumiInFiles(urls);
	goodLumiFilter.DumpToJson(((outUrl.ReplaceAll(".root",""))+".json").Data());
	}

}

