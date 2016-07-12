//
// Oleksii Toldaiev, <oleksii.toldaiev@gmail.com>
//
// ttbar to leptons event selection
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
#include "UserCode/llvv_fwk/interface/BTagCalibrationStandalone.h"

// TODO: try cmssw Btag
//#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
//#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"
// this one is for 80X -> #include "CondTools/BTau/interface/BTagCalibrationReader.h"

#include "UserCode/llvv_fwk/interface/GammaWeightsHandler.h"

#include "EgammaAnalysis/ElectronTools/interface/ElectronEnergyCalibratorRun2.h"
#include "EgammaAnalysis/ElectronTools/interface/PhotonEnergyCalibratorRun2.h" 
//#include "EgammaAnalysis/ElectronTools/interface/ElectronEnergyCalibrator.h"  
//#include "EgammaAnalysis/ElectronTools/interface/PhotonEnergyCalibrator.h" 

#include "UserCode/llvv_fwk/interface/PatUtils.h"
#include "UserCode/llvv_fwk/interface/rochcor2015.h"

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

#include <map>
#include <string>

using namespace std;

namespace utils
	{
	namespace cmssw
		{
		// TODO: it is the same jetCorrector as in MacroUtils, only Fall_ prefix is set
		// Fall15_25nsV2_
		FactorizedJetCorrector* getJetCorrector(TString baseDir, bool isMC)
			{
			gSystem->ExpandPathName(baseDir);
			//TString pf(isMC ? "MC" : "DATA");
			TString pf("Fall15_25nsV2_");
			pf += (isMC ? "MC" : "DATA");
			
			//order matters: L1 -> L2 -> L3 (-> Residuals)
			std::vector<std::string> jetCorFiles;
			std::cout << baseDir+"/"+pf+"_L1FastJet_AK4PFchs.txt" << std::endl;
			jetCorFiles.push_back((baseDir+"/"+pf+"_L1FastJet_AK4PFchs.txt").Data());
			jetCorFiles.push_back((baseDir+"/"+pf+"_L2Relative_AK4PFchs.txt").Data());
			jetCorFiles.push_back((baseDir+"/"+pf+"_L3Absolute_AK4PFchs.txt").Data());
			if(!isMC) jetCorFiles.push_back((baseDir+"/"+pf+"_L2L3Residual_AK4PFchs.txt").Data());
			// now there is a practically empty file Fall15_25nsV2_MC_L2L3Residual_AK4PFchs.txt
			// adding the run on it anyway
			//jetCorFiles.push_back((baseDir+"/"+pf+"_L2L3Residual_AK4PFchs.txt").Data());
			// it is dummy/empty file for MC and apparently is is not used
			// but in v13.1 it seemed to influence selection a bit
			// adding it for v13.4 -- will test later without it
			// and removing in 13.7 test -- compare with 13.4 & 13.4_repeat

			//init the parameters for correction
			std::vector<JetCorrectorParameters> corSteps;
			for(size_t i=0; i<jetCorFiles.size(); i++) corSteps.push_back(JetCorrectorParameters(jetCorFiles[i]));
			
			//return the corrector
			return new FactorizedJetCorrector(corSteps);
			}

		std::vector<double> smearJER(double pt, double eta, double genPt)
			{
			std::vector<double> toReturn(3,pt);
			if(genPt<=0) return toReturn;

			// These are from MacroUtils
			// FIXME: These are the 8 TeV values.
			//
			// eta=fabs(eta);
			// double ptSF(1.0), ptSF_err(0.06);
			// if(eta<0.8)                  { ptSF=1.061; ptSF_err=sqrt(pow(0.012,2)+pow(0.023,2)); }
			// else if(eta>=0.8 && eta<1.3) { ptSF=1.088; ptSF_err=sqrt(pow(0.012,2)+pow(0.029,2)); }
			// else if(eta>=1.3 && eta<1.9) { ptSF=1.106; ptSF_err=sqrt(pow(0.017,2)+pow(0.030,2)); }
			// else if(eta>=1.9 && eta<2.5) { ptSF=1.126; ptSF_err=sqrt(pow(0.035,2)+pow(0.094,2)); }
			// else if(eta>=2.5 && eta<3.0) { ptSF=1.343; ptSF_err=sqrt(pow(0.127,2)+pow(0.123,2)); }
			// else if(eta>=3.0 && eta<3.2) { ptSF=1.303; ptSF_err=sqrt(pow(0.127,2)+pow(1.303,2)); }
			// else if(eta>=3.2 && eta<5.0) { ptSF=1.320; ptSF_err=sqrt(pow(0.127,2)+pow(1.320,2)); }

			// TODO: 13TeV table for CMSSW_76X
			// https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#Smearing_procedures
			// following https://github.com/pfs/TopLJets2015/blob/master/TopAnalysis/src/CommonTools.cc
			eta=fabs(eta);
			double ptSF(1.0), ptSF_err(0.06);
			if      (eta<0.5) { ptSF=1.095; ptSF_err=0.018; }
			else if (eta<0.8) { ptSF=1.120; ptSF_err=0.028; }
			else if (eta<1.1) { ptSF=1.097; ptSF_err=0.017; }
			else if (eta<1.3) { ptSF=1.103; ptSF_err=0.033; }
			else if (eta<1.7) { ptSF=1.118; ptSF_err=0.014; }
			else if (eta<1.9) { ptSF=1.100; ptSF_err=0.033; }
			else if (eta<2.1) { ptSF=1.162; ptSF_err=0.044; }
			else if (eta<2.3) { ptSF=1.160; ptSF_err=0.048; }
			else if (eta<2.5) { ptSF=1.161; ptSF_err=0.060; }
			else if (eta<2.8) { ptSF=1.209; ptSF_err=0.059; }
			else if (eta<3.0) { ptSF=1.564; ptSF_err=0.321; }
			else if (eta<3.2) { ptSF=1.384; ptSF_err=0.033; }
			else if (eta<5.0) { ptSF=1.216; ptSF_err=0.050; }

			toReturn[0]=TMath::Max(0., (genPt+ptSF*(pt-genPt))/pt );
			toReturn[1]=TMath::Max(0., (genPt+(ptSF+ptSF_err)*(pt-genPt))/pt );
			toReturn[2]=TMath::Max(0., (genPt+(ptSF-ptSF_err)*(pt-genPt))/pt );
			return toReturn;
			}


		/* Using new stuff above
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
		//std::vector<double> smearJES(double pt, double eta, JetCorrectionUncertainty *jecUnc)
		std::vector<float> smearJES(double pt, double eta, JetCorrectionUncertainty *jecUnc)
			{
			jecUnc->setJetEta(eta);
			jecUnc->setJetPt(pt);
			double relShift=fabs(jecUnc->getUncertainty(true));
			//std::vector<double> toRet;
			std::vector<float> toRet;
			toRet.push_back((1.0+relShift)*pt);
			toRet.push_back((1.0-relShift)*pt);
			return toRet;
			}
		*/

		/* using the one in src/MacroUtils.cc
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
		*/

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
	// Set of cuts from the POG group: https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_data

	bool passID(false); 

	// float rawJetEn(jet.correctedJet("Uncorrected").energy() );

	double eta=jet.eta();
 
 	// from the twiki:
	float NHF = jet.neutralHadronEnergyFraction();
	float NEMF = jet.neutralEmEnergyFraction();
	float CHF = jet.chargedHadronEnergyFraction();
	float MUF = jet.muonEnergyFraction();
	float CEMF = jet.chargedEmEnergyFraction();
	float NumConst = jet.chargedMultiplicity()+jet.neutralMultiplicity();
	float NumNeutralParticles =jet.neutralMultiplicity();
	float CHM = jet.chargedMultiplicity(); 
	// TODO: check if these change corresponding to jet corrections? Or apply corrections after passing the selection?

	// float nhf( (jet.neutralHadronEnergy() + jet.HFHadronEnergy())/rawJetEn );
	// float nef( jet.neutralEmEnergy()/rawJetEn );
	// float cef( jet.chargedEmEnergy()/rawJetEn );
	// float chf( jet.chargedHadronEnergy()/rawJetEn );
	// float nch    = jet.chargedMultiplicity();
	// float nconst = jet.chargedMultiplicity()+jet.neutralMultiplicity();
	// float muf(jet.muonEnergy()/rawJetEn); 

	if (label=="Loose")
		// passID = ( ((nhf<0.99 && nef<0.99 && nconst>1) && ((abs(eta)<=2.4 && chf>0 && nch>0 && cef<0.99) || abs(eta)>2.4)) && abs(eta)<=3.0 );
		// the same, just with new names:
		passID = (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4) && abs(eta)<=3.0 ;
	if (label=="Tight")
		// passID = ( ((nhf<0.90 && nef<0.90 && nconst>1) && ((abs(eta)<=2.4 && chf>0 && nch>0 && cef<0.90) || abs(eta)>2.4)) && abs(eta) <=3.0);
		passID = (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4) && abs(eta)<=3.0 ;

	// there is also:
	//tightLepVetoJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || abs(eta)>2.4) && abs(eta)<=3.0

	// And the abs(eta)>3.0 part, but we never consider such jets, so... Meh!

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


/*
 * string find_W_decay(const reco::Candidate * W) {
 *
 * returns a substring {el, mu, tau, q} identifing which decay happend
 * or returns "" string if something weird happend
 *
 * Usage:
 * const reco::Candidate * W = p.daughter( W_num );
 * mc_decay += find_W_decay(W);
*/
string find_W_decay(const reco::Candidate * W) {
	const reco::Candidate * p = W; // eeh.. how C passes const arguments? are they local to function
	int d0_id, d1_id; // ids of decay daughters
	//bool found_decay=false;
	// follow W, whatever happens to it (who knows! defensive programming here)
	// until leptonic/quarkonic decay is found
	// 
	// assume there can be only W->W transitions inbetween (TODO: to check actually)
	//while (!found_decay) {
	while (true) {
		int n = p->numberOfDaughters();
		switch(n) {
		case 0: return string("");
		case 1: // it should be another W->W transition
			p = p->daughter(0);
			break; // there should not be no infinite loop here!
		case 2: // so, it should be the decay
			//found_decay = true;
			d0_id = fabs(p->daughter(0)->pdgId());
			d1_id = fabs(p->daughter(1)->pdgId());
			if (d0_id == 15 || d1_id == 15 ) return string("tau");
			if (d0_id == 11 || d1_id == 11 ) return string("el");
			if (d0_id == 13 || d1_id == 13 ) return string("mu");
			return string("q"); // FiXME: quite dangerous control-flow!
		default: // and this is just crazy
			return string("");
		}
	}
}

// int n = p.numberOfDaughters();
// if (n == 2) { // it is a decay vertes of t to something
// int d0_id = p.daughter(0)->pdgId();
// int d1_id = p.daughter(1)->pdgId();






#define MULTISEL_SIZE 256
#define NPARTICLES_SIZE 50
#define NORMKINO_DISTR_SIZE 400


// JobDef job_def = {string(isMC ? "MC": "Data"), dtag_s, job_num};

struct JobDef
{
	string type;
	string dtag;
	string job_num;
};


// inline CONTROL:

// map for TH1D and map for double (weight counters)

// inline control functions usage:
//   fill_pt_e( "control_point_name", value, weight)
//   fill_eta( "control_point_name", value, weight )   <-- different TH1F range and binning
//   increment( "control_point_name", weight )
//   printout_distrs(FILE * out)
//   printout_counters(FILE * out)

// TODO: and the multiselect weight-flow?

//std::map<string, TH1D*> th1d_distr_control;

string mc_decay("");
// Some MC datasets are inclusive, but analysis needs info on separate channels from them
// thus the events are traversed and the channel is found
// currently it is done only for TTbar channel (isTTbarMC)
//
// the sub-channel of MC is paired together with the distr_name string into the key of control distrs
// it is then printed out with dtag


std::map<std::pair <string,string>, TH1D> th1d_distr_control;
std::map<string, TH1D> th1d_distr_control_headers;
std::map<std::pair <string,string>, TH1I> th1i_distr_control;
std::map<string, TH1I> th1i_distr_control_headers;
std::map<std::pair <string,string>, double> weight_flow_control;

//std::map<string, TH1D> th1d_distr_control;
//std::map<string, TH1I> th1i_distr_control;
//std::map<string, double> weight_flow_control;


int fill_n(string control_point_name, unsigned int value, double weight)
	{
	// check if the key (mc_decay, control point) has been initialized
	std::pair <string,string> key (mc_decay, control_point_name);

	if (th1d_distr_control.find(key) == th1d_distr_control.end() )
		{
		// the control point distr has not been created/initialized
		// create it:
		//th1i_distr_control[control_point_name] = (TH1D*) new TH1D(control_point_name.c_str(), ";;Pt/E(GeV)", 400, 0., 200.);
		//th1i_distr_control[key] = TH1I(control_point_name.c_str(), ";;N", 100, 0., 100.);
		// particle counters are broken here
		// trying TH1D for v13.5
		th1d_distr_control[key] = TH1D(control_point_name.c_str(), ";;N", 100, 0., 100.);
		//cout << "creating " << control_point_name << endl;
		}

	// fill the distribution:
	// th1i_distr_control[key].Fill(value, weight);
	th1d_distr_control[key].Fill(value, weight);
	//cout << "filled " << control_point_name << endl;
	//cout << th1i_distr_control[control_point_name].Integral() << endl;

	//if (th1i_distr_control_headers.find(string("n")) == th1i_distr_control_headers.end() )
	//	{
	//	th1i_distr_control_headers[string("n")] = TH1I("Header of particle counter distributions", ";;N", 100, 0., 100.);
	//	}
	if (th1d_distr_control_headers.find(string("n")) == th1d_distr_control_headers.end() )
		{
		th1d_distr_control_headers[string("n")] = TH1D("Header of particle counter distributions", ";;N", 100, 0., 100.);
		}

	// return success:
	return 0;
	}


int fill_particle_ids(string control_point_name, int value, double weight)
	{
	// for tau (and other) fake-rates
	// check if the key (mc_decay, control point) has been initialized
	std::pair <string,string> key (mc_decay, control_point_name);

	if (th1d_distr_control.find(key) == th1d_distr_control.end() )
		{
		// the control point distr has not been created/initialized
		// create it:
		//th1i_distr_control[control_point_name] = (TH1D*) new TH1D(control_point_name.c_str(), ";;Pt/E(GeV)", 400, 0., 200.);
		//th1i_distr_control[key] = TH1I(control_point_name.c_str(), ";;N", 100, 0., 100.);
		// particle counters are broken here
		// trying TH1D for v13.5
		th1d_distr_control[key] = TH1D(control_point_name.c_str(), ";;ID", 600, -300., 300.);
		//cout << "creating " << control_point_name << endl;
		}

	// fill the distribution:
	// th1i_distr_control[key].Fill(value, weight);
	th1d_distr_control[key].Fill(value, weight);
	//cout << "filled " << control_point_name << endl;
	//cout << th1i_distr_control[control_point_name].Integral() << endl;

	if (th1d_distr_control_headers.find(string("p_id")) == th1d_distr_control_headers.end() )
		{
		th1d_distr_control_headers[string("p_id")] = TH1D("Header of particle ID distributions", ";;ID", 600, -300., 300.);
		}

	// return success:
	return 0;
	}


int fill_pu(string control_point_name, double value, double weight)
	{
	// check if the key (mc_decay, control point) has been initialized
	std::pair <string,string> key (mc_decay, control_point_name);

	if (th1d_distr_control.find(key) == th1d_distr_control.end() )
		{
		// the control point distr has not been created/initialized
		// create it:
		//th1d_distr_control[control_point_name] = (TH1D*) new TH1D(control_point_name.c_str(), ";;Pt/E(GeV)", 400, 0., 200.);
		th1d_distr_control[key] = TH1D(control_point_name.c_str(), ";;nVtx", 100, 0., 100.);
		//cout << "creating " << control_point_name << endl;
		}

	// fill the distribution:
	th1d_distr_control[key].Fill(value, weight);
	//cout << "filled " << control_point_name << endl;
	//cout << th1d_distr_control[control_point_name].Integral() << endl;

	if (th1d_distr_control_headers.find(string("pu")) == th1d_distr_control_headers.end() )
		{
		th1d_distr_control_headers[string("pu")] = TH1D("Header of pile up distributions", ";;nVtx", 100, 0., 100.);
		}

	// return success:
	return 0;
	}


int fill_pt_e(string control_point_name, double value, double weight)
	{
	// check if control point has been initialized
	std::pair <string,string> key (mc_decay, control_point_name);

	if (th1d_distr_control.find(key) == th1d_distr_control.end() )
		{
		// the control point distr has not been created/initialized
		// create it:
		//th1d_distr_control[control_point_name] = (TH1D*) new TH1D(control_point_name.c_str(), ";;Pt/E(GeV)", 400, 0., 200.);
		th1d_distr_control[key] = TH1D(control_point_name.c_str(), ";;Pt/E(GeV)", 400, 0., 400.);
		//cout << "creating " << control_point_name << endl;
		}

	// fill the distribution:
	th1d_distr_control[key].Fill(value, weight);
	//cout << "filled " << control_point_name << endl;
	//cout << th1d_distr_control[control_point_name].Integral() << endl;

	if (th1d_distr_control_headers.find(string("pt_e")) == th1d_distr_control_headers.end() )
		{
		th1d_distr_control_headers[string("pt_e")] = TH1D("Header of Pt/E distributions", ";;Pt/E(GeV)", 400, 0., 400.);
		}

	// return success:
	return 0;
	}


int fill_eta(string control_point_name, double value, double weight)
	{
	// check if control point has been initialized
	std::pair <string,string> key (mc_decay, control_point_name);

	if (th1d_distr_control.find(key) == th1d_distr_control.end() )
		{
		// the control point distr has not been created/initialized
		// create it:
		//th1d_distr_control[control_point_name] = (TH1D*) new TH1D(control_point_name.c_str(), ";;Eta", 200, -4., 4.);
		th1d_distr_control[key] = TH1D(control_point_name.c_str(), ";;Eta", 200, -4., 4.);
		}

	// fill the distribution:
	th1d_distr_control[key].Fill(value, weight);

	if (th1d_distr_control_headers.find(string("eta")) == th1d_distr_control_headers.end() )
		{
		th1d_distr_control_headers[string("eta")] = TH1D("Header of Eta distributions", ";;Eta", 200, -4., 4.);
		}

	// return success:
	return 0;
	// TODO: return the return value of Fill call?
	}




int increment(string control_point_name, double weight)
	{
	// check if control point has been initialized
	std::pair <string,string> key (mc_decay, control_point_name);

	if (weight_flow_control.find(key) == weight_flow_control.end() )
		{
		// the control point distr has not been created/initialized
		// create it:
		weight_flow_control[key] = 0.;
		}

	// fill the distribution:
	weight_flow_control[key] += weight;

	// return success:
	return 0;
	}




// Plain text output printing:

int printout_distrs(FILE * out, JobDef JD)
	{
	//th1d_distr_control
	fprintf(out, "distr_header:distr_name:header/content,data_type,dtag,job_num,distr_values\n");
	//th1d_distr_control_headers
	//th1i_distr_control_headers

	for(std::map<string, TH1D>::iterator it = th1d_distr_control_headers.begin(); it != th1d_distr_control_headers.end(); ++it)
		{
		string name = it->first;
		TH1D * header_distr = & it->second;
		// Header:
		fprintf(out, "header,%s", name.c_str());
		for (int i=0; i < header_distr->GetSize(); i++) fprintf(out, ",%g", header_distr->GetBinCenter(i));
		fprintf(out, "\n");
		}

	for(std::map<string, TH1I>::iterator it = th1i_distr_control_headers.begin(); it != th1i_distr_control_headers.end(); ++it)
		{
		string name = it->first;
		TH1I * header_distr = & it->second;
		// Header:
		fprintf(out, "header,%s", name.c_str());
		for (int i=0; i < header_distr->GetSize(); i++) fprintf(out, ",%g", header_distr->GetBinCenter(i));
		fprintf(out, "\n");
		}

	for(std::map<std::pair <string,string>, TH1D>::iterator it = th1d_distr_control.begin(); it != th1d_distr_control.end(); ++it)
		{
		// iterator->first = key
		// iterator->second = value

		const std::pair <string,string> *key = &it->first;
		string mc_decay_suffix = key->first;
		string name = key->second;

		TH1D * distr = & it->second;

		// // Content:
		fprintf(out, "%s:content, %s,%s,%s", name.c_str(), JD.type.c_str(), (JD.dtag + mc_decay_suffix).c_str(), JD.job_num.c_str());
		for (int i=0; i < distr->GetSize(); i++) fprintf(out, ",%g", distr->GetBinContent(i));
		fprintf(out, "\n");

		// New output:
		// fprintf(out, "%s:content", name.c_str());

		// JD.type.c_str(), JD.dtag.c_str(), JD.job_num.c_str()
		//for (int i=0; i < distr->GetSize(); i++)
		//	{
		//	fprintf(out, "%s,%s,%g,%g\n", name.c_str(), prefix.c_str(), distr->GetBinCenter(i), distr->GetBinContent(i));
		//	}
		//fprintf(out, "\n");
		}

	for(std::map<std::pair <string,string>, TH1I>::iterator it = th1i_distr_control.begin(); it != th1i_distr_control.end(); ++it)
		{
		// iterator->first = key
		// iterator->second = value

		const std::pair <string,string> *key = &it->first;
		string mc_decay_suffix = key->first;
		string name = key->second;

		TH1I * distr = & it->second;

		// // Content:
		fprintf(out, "%s:content, %s,%s,%s", name.c_str(), JD.type.c_str(), (JD.dtag + mc_decay_suffix).c_str(), JD.job_num.c_str());
		for (int i=0; i < distr->GetSize(); i++) fprintf(out, ",%g", distr->GetBinContent(i));
		fprintf(out, "\n");

		// New output:
		// fprintf(out, "%s:content", name.c_str());

		// JD.type.c_str(), JD.dtag.c_str(), JD.job_num.c_str()
		//for (int i=0; i < distr->GetSize(); i++)
		//	{
		//	fprintf(out, "%s,%s,%g,%g\n", name.c_str(), prefix.c_str(), distr->GetBinCenter(i), distr->GetBinContent(i));
		//	}
		//fprintf(out, "\n");
		}

	return 0;
	}


int printout_counters(FILE * out, JobDef JD)
	{
	// weight flow control
	//struct JobDef
	//{
	//	string type;
	//	string dtag;
	//	string job_num;
	//};

	// fprintf(out, "weight_flow_header:%s,%s,%s\n", JD.type.c_str(), JD.dtag.c_str(), JD.job_num.c_str());
	fprintf(out, "counter_header:val_name,data_type,dtag,job_num,value\n");

	for(std::map<std::pair <string,string>, double>::iterator it = weight_flow_control.begin(); it != weight_flow_control.end(); ++it)
		{
		const std::pair <string,string> *key = &it->first;
		string mc_decay_suffix = key->first;
		string name = key->second;

		double weight_sum = it->second;
		fprintf(out, "%s,%s,%s,%s,%g\n", name.c_str(), JD.type.c_str(), (JD.dtag + mc_decay_suffix).c_str(), JD.job_num.c_str(), weight_sum);
		}
	fprintf(out, "weight_flow end\n");

	return 0;
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

// random numbers for corrections & uncertainties
TRandom *r3 = new TRandom3();

// configure the process
const edm::ParameterSet & runProcess = edm::readPSetsFrom (argv[1])->getParameter < edm::ParameterSet > ("runProcess");

bool debug           = runProcess.getParameter<bool>  ("debug");
bool runSystematics  = runProcess.getParameter<bool>  ("runSystematics");
bool saveSummaryTree = runProcess.getParameter<bool>  ("saveSummaryTree");
bool isMC            = runProcess.getParameter<bool>  ("isMC");
double xsec          = runProcess.getParameter<double>("xsec");
int mctruthmode      = runProcess.getParameter<int>   ("mctruthmode");
TString dtag         = runProcess.getParameter<std::string>("dtag");
string dtag_s        = runProcess.getParameter<std::string>("dtag");
string job_num       = runProcess.getParameter<std::string>("job_num");

JobDef job_def = {string(isMC ? "MC": "Data"), dtag_s, job_num};

	
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

// for new orthogonality [TESTING]
bool isSingleElectronDataset = !isMC && dtag.Contains ("SingleEle");
// old trigger dataset orthogonality:
bool
	filterOnlySINGLEE  (false),
	filterOnlySINGLEMU (false);
if (!isMC)
	{
	if (dtag.Contains ("SingleMuon")) filterOnlySINGLEMU = true;
	if (dtag.Contains ("SingleEle"))  filterOnlySINGLEE  = true;
	}

// it is not used now
bool isV0JetsMC   (isMC && (dtag.Contains ("DYJetsToLL") || dtag.Contains ("WJets")));
// Reactivate for diboson shapes  
// bool isMC_ZZ      (isMC && (string (dtag.Data ()).find ("TeV_ZZ") != string::npos));
// bool isMC_WZ      (isMC && (string (dtag.Data ()).find ("TeV_WZ") != string::npos));

// adding W1Jets, W2Jets, W3Jets, W4Jets
// thus, the events with 1, 2, 3, 4 of WJets whould be excluded
// (what are the x-sections for each of these?)
bool isW0JetsSet   (isMC && (dtag.Contains ("W0Jets")));
// so, WJets provides only W0Jets events
// all the rest is in WNJets (W4Jets = W>=4Jets)
// the cross-section should probably be = WJets - sum(WNJets)
// and how to do it:
// < LHEEventProduct > lheEPHandle; lheEPHandle->hepeup().NUP  <-- this parameter
// it = 6 for W1Jets, 7 for W2Jets, 8 for W3Jets, 9 for W4Jets (why not >=9?)
// WJets have NUP = all those numbers
// for W0Jets one takes NUP = 5

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
//TString pf(isMC ? "MC" : "DATA");
//JetCorrectionUncertainty *totalJESUnc = new JetCorrectionUncertainty ((jecDir + "/"+pf+"_Uncertainty_AK4PFchs.txt").Data ());
//JetCorrectionUncertainty *totalJESUnc = new JetCorrectionUncertainty ((jecDir + "/" + (isMC ? "MC" : "DATA") + "_Uncertainty_AK4PFchs.txt").Data ());
JetCorrectionUncertainty *totalJESUnc = new JetCorrectionUncertainty ((jecDir + "/Fall15_25nsV2_" + (isMC ? "MC" : "DATA") + "_Uncertainty_AK4PFchs.txt").Data ());
// JetCorrectionUncertainty *totalJESUnc = new JetCorrectionUncertainty ((jecDir + "/MC_Uncertainty_AK4PFchs.txt").Data ());

// ------------------------------------- muon energy scale and uncertainties
// MuScleFitCorrector *muCor = NULL;
// FIXME: MuScle fit corrections for 13 TeV not available yet (more Zs are needed) getMuonCorrector (jecDir, dtag);
// TString muscleDir = runProcess.getParameter<std::string>("muscleDir");
// gSystem->ExpandPathName(muscleDir);
rochcor2015* muCor = new rochcor2015(); // This replaces the old MusScleFitCorrector that was used at RunI
// comes from:
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/RochcorMuon
// last muon POG talk:
// https://indico.cern.ch/event/533054/contributions/2171540/attachments/1274536/1891597/rochcor_run2_MuonPOG_051716.pdf

// Electron energy scale, based on https://twiki.cern.ch/twiki/bin/viewauth/CMS/EGMSmearer and adapted to this framework
string EGammaEnergyCorrectionFile = "EgammaAnalysis/ElectronTools/data/76X_16DecRereco_2015";
EpCombinationTool theEpCombinationTool;
theEpCombinationTool.init((string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/weights/GBRForest_data_25ns.root").c_str(), "gedelectron_p4combination_25ns");  //got confirmation from Matteo Sani that this works for both data and MC 
ElectronEnergyCalibratorRun2 ElectronEnCorrector(theEpCombinationTool, isMC, false, EGammaEnergyCorrectionFile);
ElectronEnCorrector.initPrivateRng(new TRandom(1234));


// --------------------------------------- lepton efficiencies
LeptonEfficiencySF lepEff;

// --------------------------------------- b-tagging 
// TODO: move all these numbers to where they are applied??
// btagMedium is used twice in the code
// merge those tagging procedures and eliminated the variable?

// Prescriptions taken from: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation74X
// TODO: update to https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation76X

// b-tagging working points for 50ns 
//   (pfC|c)ombinedInclusiveSecondaryVertexV2BJetTags
//      v2CSVv2L 0.605
//      v2CSVv2M 0.890
//      v2CSVv2T 0.970
double
	btagLoose(0.605), // not used anywhere in the code
	//btagMedium(0.890), // used twice in the code
	btagMedium(0.8), // new medium working point
	btagTight(0.970); // not used anywhere in the code

//b-tagging: scale factors
//beff and leff must be derived from the MC sample using the discriminator vs flavor
//the scale factors are taken as average numbers from the pT dependent curves see:
//https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG#2012_Data_and_MC_EPS13_prescript
BTagSFUtil btsfutil;
double beff(0.68), sfb(0.99), sfbunc(0.015);
double leff(0.13), sfl(1.05), sflunc(0.12);

// Btag SF and eff from https://indico.cern.ch/event/437675/#preview:1629681
sfb = 0.861; // SF is not used --- BTagCalibrationReader btagCal instead
// sbbunc =;
beff = 0.559;

// new btag calibration
// TODO: check callibration readers in https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration
// and latest standalone callibrator:
// https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/CondTools/BTau/test/BTagCalibrationStandalone.h

// Setup calibration readers
//BTagCalibration btagCalib("CSVv2", string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/weights/btagSF_CSVv2.csv");
BTagCalibration btagCalib("CSVv2", string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/weights/CSVv2_76X.csv");
// TODO: try new llvv_fwk/data/weights/CSVv2_76X.csv
// and there:
//The name of the measurements is
//
//  * "incl" for light jets,
//  * "mujets" for b and c jets for what concerns the pT/eta dependence for the different WP for JP and CSVv2 and
//  * "ttbar" for b and c jets for what concerns the pT/eta dependence for the different WP for cMVAv2, but only to be used for jets with a pT spectrum similar to that in ttbar.
//  * The measurement "iterativefit" provides the SF as a function of the discriminator shape. 
// --- so "incl" instead of "comb" for light-quarks

// TODO: update btag CSVv2
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation76X#Data_MC_Scale_Factors
// OP_MEDIUM instead of OP_LOOSE?
BTagCalibrationReader btagCal   (&btagCalib, BTagEntry::OP_MEDIUM, "mujets", "central");  // calibration instance, operating point, measurement type, systematics type
BTagCalibrationReader btagCalUp (&btagCalib, BTagEntry::OP_MEDIUM, "mujets", "up"     );  // sys up
BTagCalibrationReader btagCalDn (&btagCalib, BTagEntry::OP_MEDIUM, "mujets", "down"   );  // sys down
//BTagCalibrationReader btagCalL  (&btagCalib, BTagEntry::OP_LOOSE, "comb", "central");  // calibration instance, operating point, measurement type, systematics type
//BTagCalibrationReader btagCalLUp(&btagCalib, BTagEntry::OP_LOOSE, "comb", "up"     );  // sys up
//BTagCalibrationReader btagCalLDn(&btagCalib, BTagEntry::OP_LOOSE, "comb", "down"   );  // sys down
BTagCalibrationReader btagCalL  (&btagCalib, BTagEntry::OP_MEDIUM, "incl", "central");  // calibration instance, operating point, measurement type, systematics type
BTagCalibrationReader btagCalLUp(&btagCalib, BTagEntry::OP_MEDIUM, "incl", "up"     );  // sys up
BTagCalibrationReader btagCalLDn(&btagCalib, BTagEntry::OP_MEDIUM, "incl", "down"   );  // sys down



/* TODO: CMSSW calibration:
The twiki https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation76X#Supported_Algorithms_and_Operati says

The name of the measurements is
 * "incl" for light jets,
 * "mujets" for b and c jets for what concerns the pT/eta dependence for the different WP for JP and CSVv2 and 
 * "ttbar" for b and c jets for what concerns the pT/eta dependence for the different WP for cMVAv2, but only to be used for jets with a pT spectrum similar to that in ttbar.
 * The measurement "iterativefit" provides the SF as a function of the discriminator shape. 


// setup calibration readers
// Before CMSSW_8_0_11 and CMSSW_8_1_0, one reader was needed 
// for every OperatingPoint/MeasurementType/SysType combination.
BTagCalibration calib("csvv1", "CSVV1.csv");
BTagCalibrationReader reader(BTagEntry::OP_LOOSE,  // operating point
                             "central");           // systematics type
BTagCalibrationReader reader_up(BTagEntry::OP_LOOSE, "up");  // sys up
BTagCalibrationReader reader_do(BTagEntry::OP_LOOSE, "down");  // sys down

reader.load(&calib,               // calibration instance
            BTagEntry::FLAV_B,    // btag flavour
            "comb")               // measurement type
// reader_up.load(...)
// reader_down.load(...)



// Usage:

float JetPt = b_jet.pt(); bool DoubleUncertainty = false;
if (JetPt>MaxBJetPt)  { // use MaxLJetPt for  light jets
  JetPt = MaxBJetPt; 
  DoubleUncertainty = true;
}  

// Note: this is for b jets, for c jets (light jets) use FLAV_C (FLAV_UDSG)
double jet_scalefactor = reader.eval(BTagEntry::FLAV_B, b_jet.eta(), JetPt); 
double jet_scalefactor_up =  reader_up.eval(BTagEntry::FLAV_B, b_jet.eta(), JetPt); 
double jet_scalefactor_do =  reader_do.eval(BTagEntry::FLAV_B, b_jet.eta(), JetPt); 

if (DoubleUncertainty) {
   jet_scalefactor_up = 2*(jet_scalefactor_up - jet_scalefactor) + jet_scalefactor; 
   jet_scalefactor_do = 2*(jet_scalefactor_do - jet_scalefactor) + jet_scalefactor; 
}

*/




// ------------------------------ electron IDs
// does not appear anywhere in the code at all
// TString
	// electronIdMainTag("cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
	// electronIdVetoTag("cutBasedElectronID-Spring15-25ns-V1-standalone-tight");

// --------------------------------------- pileup weighting
// pile-up is done directly with direct_pileup_reweight
std::vector<double> direct_pileup_reweight = runProcess.getParameter < std::vector < double >>("pileup_reweight_direct");
std::vector<double> direct_pileup_reweight_up = runProcess.getParameter < std::vector < double >>("pileup_reweight_direct_up");
std::vector<double> direct_pileup_reweight_down = runProcess.getParameter < std::vector < double >>("pileup_reweight_direct_down");
	
gROOT->cd ();                 //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE

// --------------------------------------- hardcoded MET filter
// trying met-filters in HLT paths (76X MINIAODs should have it)
/*
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
*/

// ----------------------------
// So here we got all the parameters from the config




//##############################################
//########    INITIATING HISTOGRAMS     ########
//##############################################

// Removed the SmartSelectionMonitor
// SmartSelectionMonitor mon;

TH1D* singlelep_ttbar_initialevents  = (TH1D*) new TH1D("singlelep_ttbar_init",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  ); 
TH1D* singlelep_ttbar_preselectedevents = (TH1D*) new TH1D("singlelep_ttbar_presele",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  ); 
/* These counter histograms are dissabled -- use int for that
TH1D* singlelep_ttbar_selected_mu_events = (TH1D*) new TH1D("singlelep_ttbar_sele_mu",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  ); 
TH1D* singlelep_ttbar_selected_el_events = (TH1D*) new TH1D("singlelep_ttbar_sele_el",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  ); 
TH1D* singlelep_ttbar_selected2_mu_events = (TH1D*) new TH1D("singlelep_ttbar_sele2_mu",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  ); 
TH1D* singlelep_ttbar_selected2_el_events = (TH1D*) new TH1D("singlelep_ttbar_sele2_el",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  ); 

TH1D* singlelep_ttbar_maraselected_mu_events = (TH1D*) new TH1D("singlelep_ttbar_marasele_mu",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  ); 
TH1D* singlelep_ttbar_maraselected_el_events = (TH1D*) new TH1D("singlelep_ttbar_marasele_el",     ";Transverse momentum [GeV];Events",            100, 0.,  500.  ); 
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
string FileName = ((outUrl.ReplaceAll(".root",""))+".csv").Data();
csv_out = fopen(FileName.c_str(), "w");

// FIXME: it's initialization of a rare control point, make it automatic somehow? initialize all?
// to check if processing goes well now
// increment( string("weight_passed_oursel"), 0. );

fprintf(csv_out, "Headers\n");

// MULTISELECT headers:

fprintf(csv_out, "weights_in_no_channel, num_events, num_events_pass_lumi, initial_events_weight, sum_rawweights_pass_trig, sum_weights_passtrig, weight_before_channel_select");
for (int i=0; i<MULTISEL_SIZE; i++)
	{
	fprintf(csv_out, ",%d", i);
	}
fprintf(csv_out, "\n");

fprintf(csv_out, "weights_in_el_channel, num_events, num_events_pass_lumi, initial_events_weight, sum_rawweights_pass_trig, sum_weights_passtrig, weight_before_channel_select");
for (int i=0; i<MULTISEL_SIZE; i++)
	{
	fprintf(csv_out, ",%d", i);
	}
fprintf(csv_out, "\n");

fprintf(csv_out, "weights_in_mu_channel, num_events, num_events_pass_lumi, initial_events_weight, sum_rawweights_pass_trig, sum_weights_passtrig, weight_before_channel_select");
for (int i=0; i<MULTISEL_SIZE; i++)
	{
	fprintf(csv_out, ",%d", i);
	}
fprintf(csv_out, "\n");

fprintf(csv_out, "weights_in_elmu_channel, num_events, num_events_pass_lumi, initial_events_weight, sum_rawweights_pass_trig, sum_weights_passtrig, weight_before_channel_select");
for (int i=0; i<MULTISEL_SIZE; i++)
	{
	fprintf(csv_out, ",%d", i);
	}
fprintf(csv_out, "\n");

fprintf(csv_out, "weights_in_elel_channel, num_events, num_events_pass_lumi, initial_events_weight, sum_rawweights_pass_trig, sum_weights_passtrig, weight_before_channel_select");
for (int i=0; i<MULTISEL_SIZE; i++)
	{
	fprintf(csv_out, ",%d", i);
	}
fprintf(csv_out, "\n");

fprintf(csv_out, "weights_in_mumu_channel, num_events, num_events_pass_lumi, initial_events_weight, sum_rawweights_pass_trig, sum_weights_passtrig, weight_before_channel_select");
for (int i=0; i<MULTISEL_SIZE; i++)
	{
	fprintf(csv_out, ",%d", i);
	}
fprintf(csv_out, "\n");



// EVENTOUT

fprintf(csv_out, "crossel:pu_num_inters,nGoodPV, rawWeight, weight, isElectron, l_px,l_py,l_pz,l_e, b1_px,b1_py,b1_pz,b1_e, j1_px,j1_py,j1_pz,j1_e,j2_px,j2_py,j2_pz,j2_e\n");
fprintf(csv_out, "oursel: iev, pu_num_inters,nGoodPV, rawWeight, weight, weight_up, weight_down, isElectron,");
fprintf(csv_out, "met_0, met_1, met_2, met_3, met_4, met_5, met_6,");
fprintf(csv_out, "l_vz, l_px,l_py,l_pz,l_e,");
fprintf(csv_out, "tau_vz, tau_px,tau_py,tau_pz,tau_e,");
fprintf(csv_out, "b1_vz, b1_px,b1_py,b1_pz,b1_e,");
fprintf(csv_out, "j1_vz, j1_px,j1_py,j1_pz,j1_e,");
fprintf(csv_out, "j2_vz, j2_px,j2_py,j2_pz,j2_e\n");
//fprintf(csv_out, "j1_vz, j1_pt, j1_pt_up, j1_pt_down, j1_px,j1_py,j1_pz,j1_e,");
//fprintf(csv_out, "j2_vz, j2_pt, j2_pt_up, j2_pt_down, j2_px,j2_py,j2_pz,j2_e\n");

fprintf(csv_out, "marasel:pu_num_inters,nGoodPV, rawWeight, weight, isElectron, l_px,l_py,l_pz,l_e, b1_px,b1_py,b1_pz,b1_e,b2_px,b2_py,b2_pz,b2_e, j1_px,j1_py,j1_pz,j1_e,j2_px,j2_py,j2_pz,j2_e,j3_px,j3_py,j3_pz,j3_e,j4_px,j4_py,j4_pz,j4_e\n");

fprintf(csv_out, "\n");

for(size_t f=0; f<urls.size();++f)
	{
	fprintf(csv_out, "Processing file: %s\n", urls[f].c_str());
	cout << "Processing file: " << urls[f].c_str() << "\n";
	TFile* file = TFile::Open(urls[f].c_str());
	fwlite::Event ev(file);
	printf ("Scanning the ntuple %2lu/%2lu :",f+1, urls.size());

	// acceptance parameters
	unsigned int iev(0); // number of events
	unsigned int n_events_pass_lumi(0); // number of events lassing lumi

	double sum_weights_raw = 0; // sum of raw weights
	double sum_weights = 0; // sum of final weights

	// sum weights before the particle selection
	double sum_weights_passtrig_raw = 0;
	double sum_weights_passtrig = 0;

	// before channel multiselect
	double weight_before_channel_select = 0;

	double crossel_sum_weights_raw = 0; // crossel
	double crossel_sum_weights = 0;
	double oursel_sum_weights_raw = 0; // oursel
	double oursel_sum_weights = 0;
	double oursel_sum_weights_el = 0;
	double oursel_sum_weights_mu = 0;
	double marasel_sum_weights_raw = 0; // marasel
	double marasel_sum_weights = 0;



	// MULTISELECTION array -- 8 bits
	// now 6 bits are used -- 0-63 is max
	double weights_in_no_channel[MULTISEL_SIZE], weights_in_el_channel[MULTISEL_SIZE], weights_in_mu_channel[MULTISEL_SIZE],
		weights_in_elmu_channel[MULTISEL_SIZE], weights_in_elel_channel[MULTISEL_SIZE], weights_in_mumu_channel[MULTISEL_SIZE];
	//int weights_in_selections_int[100];
	for (int i=0; i<MULTISEL_SIZE; i++)
		{
		//weights_in_selections[i] = 0;
		weights_in_no_channel[i] = 0;
		weights_in_el_channel[i] = 0;
		weights_in_mu_channel[i] = 0;
		weights_in_elmu_channel[i] = 0;
		weights_in_elel_channel[i] = 0;
		weights_in_mumu_channel[i] = 0;
		//weights_in_selections_int[i] = 0;
		}

	unsigned int negative_event_nvtx[100];
	unsigned int positive_event_nvtx[100];
	double negative_event_pernvtx_weight[100];
	double positive_event_pernvtx_weight[100];
	double negative_event_pergoodpv_weight[100];
	double positive_event_pergoodpv_weight[100];
	double event_pergoodpv_weight[100];
	double n_selected_leptons_weighted[100];
	double n_selected_taus_weighted[100];
	double n_selected_jets_weighted[100];
	double n_selected_bjets_weighted[100];
	for (int i=0; i<100; i++)
		{
		negative_event_nvtx[i] = 0;
		positive_event_nvtx[i] = 0;
		negative_event_pernvtx_weight[i] = 0;
		positive_event_pernvtx_weight[i] = 0;
		negative_event_pergoodpv_weight[i] = 0;
		positive_event_pergoodpv_weight[i] = 0;
		event_pergoodpv_weight[i] = 0;
		n_selected_leptons_weighted[i] = 0;
		n_selected_taus_weighted[i] = 0;
		n_selected_jets_weighted[i] = 0;
		n_selected_bjets_weighted[i] = 0;
		}

	int treeStep (ev.size()/50);

	for (ev.toBegin(); !ev.atEnd(); ++ev)
		{
		if(debug && iev == 2){
			cout << "Finished processing the " << iev << " event in the file, exiting" << endl;
			//return 0;
			break;
			}

		mc_decay = string("");

		reco::GenParticleCollection gen;
		fwlite::Handle<reco::GenParticleCollection> genHandle;
		genHandle.getByLabel(ev, "prunedGenParticles");
		if(genHandle.isValid() ) gen = *genHandle;

		fwlite::Handle < LHEEventProduct > lheEPHandle;
		lheEPHandle.getByLabel (ev, "externalLHEProducer");

		if (debug)
			{
			cout << "number of gen particles = " << gen.size() << "\n";

			//mon.fillHisto ("nup", "", lheEPHandle->hepeup ().NUP, 1);
			cout << "nup = " << lheEPHandle->hepeup().NUP << "\n";
			}

		// take only W0Jets events from WJets set: (W0Jets have NUP == 5)
		if (isW0JetsSet && (lheEPHandle->hepeup().NUP != 5))
			continue;

		// -------------------------- trying to extract what decay was generated here
		// iEvent.getByLabel("genParticles", genParticles);
		// for(size_t i = 0; i < genParticles->size(); ++ i) {
		// }
		// listings of particles' mother-daughters  show
		// W status 44 can turn into the same W status 44 several times
		// W status 62 decays into leton-neutrino
		//
		// otozh
		// there are numerous ways to do it:
		//  - find motherless particles (protons in our case) and move from them
		//  - just check out two first particles in the collection -- it seems these are the mothers always
		//  - find t-quarks among particles and transverse their branches,
		//    hoping that t and tbar cannot happen in not a ttbar pair
		//
		// let's do the naive 3 version
		// targeting:
		//   find decay vertex t-W
		//   follow along W and find decay vertex to leptons or quarks
		// -- thus recursion is needed

		// traversing separate t quarks
		//if (debug) {
		if (isTTbarMC && genHandle.isValid()) {
			//string mc_decay(""); // move to all job parameters
			// every found t-quark decay branch will be added as a substring to this string (el, mu, tau, q)
			// TODO: what to do if there are > t-decays than 1 ttbar pair in an event? (naive traverse approach also doesn't work?)
			for(size_t i = 0; i < gen.size(); ++ i) {
				const reco::GenParticle & p = gen[i];
				int id = p.pdgId();
				int st = p.status(); // TODO: check what is status in decat simulation (pythia for our TTbar set)

				if (id == 6) { // if it is a t quark
					// looking for W+ and its' leptonic decay
					// without checks, just checking 2 daughters
					int n = p.numberOfDaughters();
					if (n == 2) { // it is a decay vertes of t to something
						int d0_id = p.daughter(0)->pdgId();
						int d1_id = p.daughter(1)->pdgId();
						int W_num = d0_id == 24 ? 0 : (d1_id == 24 ? 1 : -1) ;
						if (W_num < 0) continue;
						const reco::Candidate * W = p.daughter( W_num );
						mc_decay += find_W_decay(W);
					}
				}

				if (id == -6) { // if it is a tbar quark
					// looking for W- and its' leptonic decay
					int n = p.numberOfDaughters();
					if (n == 2) { // it is a decay vertes of t to something
						int d0_id = p.daughter(0)->pdgId();
						int d1_id = p.daughter(1)->pdgId();
						int W_num = d0_id == -24 ? 0 : (d1_id == -24 ? 1 : -1) ;
						if (W_num < 0) continue;
						const reco::Candidate * W = p.daughter( W_num );
						mc_decay += find_W_decay(W);
						mc_decay += string("bar");
					}
				}
			}
			// so, mc_decay will be populated with strings matching t decays
			// hopefully, in TTbar sample only ttbar decays are present and it is not ambigous
		}

		if (debug) {
			cout << "MC suffix " << mc_decay << " is found\n";
			}

		if (!mc_decay.empty()) mc_decay = string("_") + mc_decay;

		//* List of mother-daughters for all particles
		//* TODO: make it into a separate function

		/*
		if (debug) {
			for(size_t i = 0; i < gen.size(); ++ i) {
				const reco::GenParticle & p = gen[i];
				int id = p.pdgId();
				int st = p.status();
				int n = p.numberOfDaughters();
				cout << i << ": " << id << " " << st;
				if (p.numberOfMothers() != 0) cout << " <- " ;
				for (int j = 0 ; j < p.numberOfMothers(); ++j) {
					const reco::Candidate * mom = p.mother(j);
					cout << " " << mom->pdgId() << " " << mom->status() << ";";
					}
				cout << "\n";
				if (n>0) {
					cout << "\t|-> " ;
					for (int j = 0; j < n; ++j) {
						const reco::Candidate * d = p.daughter( j );
						cout << d->pdgId() << " " << d->status() << "; " ;
						}
					cout << "\n";
					}
				}
			}
		*/

		// inline control functions usage:
		//   fill_pt_e( "control_point_name", value, weight)
		//   fill_eta( "control_point_name", value, weight )   <-- different TH1F range and binning
		//   increment( "control_point_name", weight )

		increment( string("weightflow_n_miniaod_events"), 1.0 );

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
		double puWeight_up      (1.0);
		double puWeight_down    (1.0);

		// rawWeight is everything but Pile-Up
		double rawWeight        (1.0);

		double HLT_efficiency_sf = 1.0;

		// final weight of the event
		double weight           (1.0);
		double weight_up        (1.0);
		double weight_down      (1.0);
		// and systematic corrections? TODO: check how TotalWeight_plus is used?


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
		weight_up *= weightGen;
		weight_down *= weightGen;
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
		unsigned int num_inters = 0, num_inters_raw = 0;
		double weight_pu_test = weight;
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
			if (num_inters<100) {
				puWeight = direct_pileup_reweight[num_inters];
				puWeight_up = direct_pileup_reweight_up[num_inters];
				puWeight_down = direct_pileup_reweight_down[num_inters];
			}
			else {//puWeight = 0; puWeight_up = 0; puWeight_down = 0;
				continue; // just move on
			}

			// FIXME: testing raw vtx.size()
			num_inters_raw = vtx.size();
			if (num_inters_raw<100) {weight_pu_test *= direct_pileup_reweight[num_inters_raw];}
			else {weight_pu_test *= 1.5e-16;}
			// TODO: implement error margins of pile-up
			}
		else
			{
			// get data pile-up into num_inters
			// for now use number of good vertices for the data
			// should substitute it with something more appropriate
			// num_inters = nGoodPV;
			// let's try using the size of primary vertex collection (before selecting the good vertices)
			num_inters = vtx.size();
			}

		weight *= puWeight;
		weight_up *= puWeight_up;
		weight_down *= puWeight_down;

		// --------------- here the weighting/shaping of MC should be done
		// --------------------- save distributions of weights
		sum_weights += weight;
		sum_weights_raw += rawWeight;

		// inline control functions usage:
		//   fill_pt_e( "control_point_name", value, weight)
		//   fill_eta( "control_point_name", value, weight )   <-- different TH1F range and binning
		//   increment( "control_point_name", weight )

		increment( string("weightflow_weighted_raw_miniaod_events"), rawWeight );
		increment( string("weightflow_weighted_miniaod_events"), weight );
		increment( string("weightflow_weighted_up_miniaod_events"), weight_up );
		increment( string("weightflow_weighted_down_miniaod_events"), weight_down );

		fill_pu( string("pileup_rawweight_perrawvtxsize"), vtx.size(), rawWeight);
		fill_pu( string("pileup_weight_perrawvtxsize"), vtx.size(), weight_pu_test);

		fill_pu( string("pileup_ini_rawweight_pergoodpv"), nGoodPV, rawWeight);
		fill_pu( string("pileup_ini_weight_pergoodpv"), nGoodPV, weight);
		fill_pu( string("pileup_ini_weight_up_pergoodpv"), nGoodPV, weight_up);
		fill_pu( string("pileup_ini_weight_down_pergoodpv"), nGoodPV, weight_down);

		fill_pu( string("pileup_rawweight_pernuminters"), num_inters, rawWeight);
		fill_pu( string("pileup_weight_pernuminters"), num_inters, weight);
		fill_pu( string("pileup_weight_up_pernuminters"), num_inters, weight_up);
		fill_pu( string("pileup_weight_down_pernuminters"), num_inters, weight_down);

		//num_inters = 1;
		if (num_inters>99) num_inters = 99;
		if (nGoodPV>100) nGoodPV = 99;
		event_pergoodpv_weight[nGoodPV] += weight;
		//if (num_inters<0)  num_inters = 0;
		if (weightGen<0)
			{
			increment( string("negative_events"), 1 );
			fill_pu( string("pileup_negative_weight_pernuminters"), num_inters, weight);
			fill_pu( string("pileup_negative_rawweight_pernuminters"), num_inters, rawWeight);

			negative_event_nvtx[num_inters] += 1;
			negative_event_pernvtx_weight[num_inters] += weight;
			negative_event_pergoodpv_weight[nGoodPV] += weight;
			}
		else
			{
			increment( string("positive_events"), 1 );
			fill_pu( string("pileup_positive_weight_pernuminters"), num_inters, weight);
			fill_pu( string("pileup_positive_rawweight_pernuminters"), num_inters, rawWeight);
			positive_event_nvtx[num_inters] += 1;
			positive_event_pernvtx_weight[num_inters] += weight;
			positive_event_pergoodpv_weight[nGoodPV] += weight;
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
		if(!goodLumiFilter.isGoodLumi(ev.eventAuxiliary().run(), ev.eventAuxiliary().luminosityBlock())) continue; 
		// Notice: it is the first continue in the event loop
		n_events_pass_lumi += 1;
		// there is no sum_weights_pass_lumi -- lumi is for data only..

		// inline control functions usage:
		//   fill_pt_e( "control_point_name", value, weight)
		//   fill_eta( "control_point_name", value, weight )   <-- different TH1F range and binning
		//   increment( "control_point_name", weight )

		increment( string("weightflow_weight_passed_lumi"), weight ); // should not matter
		increment( string("weightflow_weight_up_passed_lumi"), weight_up ); // should not matter
		increment( string("weightflow_weight_down_passed_lumi"), weight_down ); // should not matter


		// --------------------------------------------- apply trigger
		// ---------------- and require compatibilitiy of the event with the PD
		edm::TriggerResultsByName tr = ev.triggerResultsByName ("HLT");
		if (!tr.isValid ()){
			cout << "Trigger is not valid" << endl;
			return false;
			}

		if(debug){
			cout << "Printing HLT trigger list" << endl;
			//cout << "-- Commented out --" << endl;
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
		
		//if(filterOnlySINGLEMU) {                    eTrigger = false; }
		//if(filterOnlySINGLEE)  { muTrigger = false;                   }
		// SingleMuon-dataset jobs process double-HLT events
		// SingleElectron-dataset jobs skip them

		// if data and SingleElectron dataset and both triggers -- skip event
		if (!debug) {
			if (!isMC && isSingleElectronDataset && eTrigger && muTrigger) continue;
		
			if (!(eTrigger || muTrigger)) continue;   //ONLY RUN ON THE EVENTS THAT PASS OUR TRIGGERS
		}

		// TODO: ----------------------------- HLT efficiency scale factors
		// https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffsRun2#Muon_reconstruction_identificati
		// -- SingleMu Triggers
		// ??? https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2#Efficiencies_and_scale_factors
		// ??? -- Triggering electrons MVA-based WPs, scale factors for 76X
		//
		// in principle, one needs the electron/muon which triggered the HLT
		// and according to its' pt, eta the scale factor is extracted
		//
		// then, how should one proceed?
		// if there are both HLTs on, should both factors be applied?

		// double HLT_efficiency_sf = 1.0;

		//HLT_efficiency_sf *= eTrigger  ? eHLT_sf[] : 1 ;
		//HLT_efficiency_sf *= muTrigger ? muHLT_SF[] : 1 ;


		sum_weights_passtrig_raw += rawWeight;
		sum_weights_passtrig += weight;


		// inline control functions usage:
		//   fill_pt_e( "control_point_name", value, weight)
		//   fill_eta( "control_point_name", value, weight )   <-- different TH1F range and binning
		//   increment( "control_point_name", weight )

		increment( string("weightflow_weight_passed_trig"), weight ); // should not matter
		increment( string("weightflow_weight_up_passed_trig"), weight_up ); // should not matter
		increment( string("weightflow_weight_down_passed_trig"), weight_down ); // should not matter

		if(debug)
			{
			cout << "Set triggers" << endl;
			}

		fill_pu( string("pileup_passtrig_rawweight_pergoodpv"), nGoodPV, rawWeight);
		fill_pu( string("pileup_passtrig_weight_pergoodpv"), nGoodPV, weight);
		fill_pu( string("pileup_passtrig_weight_up_pergoodpv"), nGoodPV, weight_up);
		fill_pu( string("pileup_passtrig_weight_down_pergoodpv"), nGoodPV, weight_down);

		fill_pu( string("pileup_passtrig_rawweight_perrawvtxsize"), vtx.size(), rawWeight);
		fill_pu( string("pileup_passtrig_weight_perrawvtxsize"), vtx.size(), weight_pu_test);
		fill_pu( string("pileup_passtrig_rawweight_pernuminters"), num_inters, rawWeight);
		fill_pu( string("pileup_passtrig_weight_pernuminters"), num_inters, weight);


		// ------------------------------------------------- Apply MET filters
		//if( !isMC && !metFiler.passMetFilter( ev, isPromptReco)) continue;
		// New passMetFilter procedure from patUtils:
		//if( !isMC && !metFiler.passMetFilter( ev )) continue;
		// trying the HLTs path for metfilters:
		// https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2015#ETmiss_filters
		// https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#MiniAOD_76X_v2_produced_with_the
		// recommendations (6-2016):
		// Flag_HBHENoiseFilter TO BE USED
		// Flag_HBHENoiseIsoFilter TO BE USED
		// Flag_CSCTightHalo2015Filter TO BE USED
		// Flag_EcalDeadCellTriggerPrimitiveFilter TO BE USED
		// Flag_goodVertices TO BE USED
		// Flag_eeBadScFilter TO BE USED 
		edm::TriggerResultsByName metFilters = ev.triggerResultsByName("PAT");   //is present only if PAT (and miniAOD) is not run simultaniously with RECO
		if(!metFilters.isValid()){metFilters = ev.triggerResultsByName("RECO");} //if not present, then it's part of RECO
		if(!metFilters.isValid()){       
			printf("TriggerResultsByName for MET filters is not found in the process, as a consequence the MET filter is disabled for this event\n");    
		}

		if (! isMC && metFilters.isValid()) {
			if(debug){
				cout << "Printing PAT/RECO trigger list for MET filters here" << endl;
				for(edm::TriggerNames::Strings::const_iterator trnames = metFilters.triggerNames().begin(); trnames!=metFilters.triggerNames().end(); ++trnames)
					cout << *trnames << endl;
				cout << "----------- End of trigger list ----------" << endl;
				//return 0;
			}
			//std::vector<std::string>& patterns("Flag_HBHENoiseFilter", "Flag_HBHENoiseIsoFilter", "Flag_CSCTightHalo2015Filter", "Flag_EcalDeadCellTriggerPrimitiveFilter", "Flag_goodVertices", "Flag_eeBadScFilter");
			if (!utils::passTriggerPatterns(metFilters, "Flag_HBHENoiseFilter*", "Flag_HBHENoiseIsoFilter*", "Flag_CSCTightHalo2015Filter*", "Flag_EcalDeadCellTriggerPrimitiveFilter*"))
				continue;
			if (!utils::passTriggerPatterns(metFilters, "Flag_goodVertices")) continue;
			if (!utils::passTriggerPatterns(metFilters, "Flag_eeBadScFilter")) continue;
		}
		

		if(debug)
			{
			cout << "met filters applied here" << endl;
			}


		// ------------------------- event physics and the corresponding selection

		//------------------------- load all the objects we will need to access

		double rho = 0;
		fwlite::Handle<double> rhoHandle;
		rhoHandle.getByLabel(ev, "fixedGridRhoFastjetAll");
		if(rhoHandle.isValid() ) rho = *rhoHandle;




		// FIXME: Save time and don't load the rest of the objects when selecting by mctruthmode :)
		//bool hasTop(false);
		//int
		//ngenLeptonsStatus3(0),
		//ngenLeptonsNonTauSonsStatus3(0),
		//ngenTausStatus3(0),
		//ngenQuarksStatus3(0);

		//double tPt(0.), tbarPt(0.); // top pt reweighting - dummy value results in weight equal to 1 if not set in loop
		//float wgtTopPt(1.0), wgtTopPtUp(1.0), wgtTopPtDown(1.0);
		// TODO: what is this??
		// there was some wague answer from Pietro.....
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
				if (
					(ngenLeptonsNonTauSonsStatus3!=1 || ngenTausStatus3!=1 ) &&
					(ngenLeptonsNonTauSonsStatus3!=2                      ) &&
					(ngenLeptonsNonTauSonsStatus3+ngenTausStatus3!=1      ) 
					)
				isHad=true;
					
				//if(mctruthmode==6 && (ngenLeptonsNonTauSonsStatus3!=0 || ngenTausStatus3!=0  || !hasTop )) continue;
				if(mctruthmode==6 && (!isHad || !hasTop )) continue;
				}
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

		pat::JetCollection jets;
		fwlite::Handle<pat::JetCollection>jetsHandle;
		jetsHandle.getByLabel(ev, "slimmedJets");
		if(jetsHandle.isValid() ) jets = *jetsHandle;

		// testing WNJets
		if(debug){
			cout << "N Jets:" << jets.size() << "\n";
			}

		pat::PhotonCollection photons;
		fwlite::Handle<pat::PhotonCollection> photonsHandle;
		photonsHandle.getByLabel(ev, "slimmedPhotons");
		if(photonsHandle.isValid() ) photons = *photonsHandle;

		pat::METCollection mets;
		fwlite::Handle<pat::METCollection> metsHandle;
		metsHandle.getByLabel(ev, "slimmedMETs");
		if(metsHandle.isValid() ) mets = *metsHandle;
		pat::MET met = mets[0];
		// LorentzVector met = mets[0].p4 ();


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


		// ------------------------------------ merging electrons and muons
		// Let's merge after processing and channel assignment
		// std::vector<patUtils::GenericLepton> leptons;
		// for(size_t l=0; l<electrons.size(); ++l) leptons.push_back(patUtils::GenericLepton (electrons[l] ));
		// for(size_t l=0; l<muons.size(); ++l)     leptons.push_back(patUtils::GenericLepton (muons[l]     ));
		// std::sort(leptons.begin(), leptons.end(), utils::sort_CandidatesByPt);


		// CONTROLINFO
		// Control values for raw particles:

		std::sort (muons.begin(),  muons.end(),  utils::sort_CandidatesByPt);
		for(size_t n=0; n<muons.size(); ++n)
			{
			fill_pt_e( string("all_muons_slimmed_pt"), muons[n].pt(), weight);
			if (n < 2)
				{
				fill_pt_e( string("top2pt_muons_slimmed_pt"), muons[n].pt(), weight);
				}
			}

		std::sort (electrons.begin(),  electrons.end(),  utils::sort_CandidatesByPt);
		for(size_t n=0; n<electrons.size(); ++n)
			{
			fill_pt_e( string("all_electrons_slimmed_pt"), electrons[n].pt(), weight);
			if (n < 2)
				{
				fill_pt_e( string("top2pt_electrons_slimmed_pt"), electrons[n].pt(), weight);
				}
			}


		std::sort (taus.begin(),  taus.end(),  utils::sort_CandidatesByPt);
		for(size_t n=0; n<taus.size(); ++n)
			{
			fill_pt_e( string("all_taus_slimmed_pt"), taus[n].pt(), weight);
			if (n < 1)
				{
				fill_pt_e( string("top1pt_taus_slimmed_pt"), taus[n].pt(), weight);
				}
			}

		std::sort (jets.begin(),  jets.end(),  utils::sort_CandidatesByPt);
		for(size_t n=0; n<jets.size(); n++)
		{
			fill_pt_e( string("all_jets_slimmed_pt"), jets[n].pt(), weight);
			if (n < 2)
				{
				fill_pt_e( string("top2pt_jets_slimmed_pt"), jets[n].pt(), weight);
				}
		}



		// and also mets:
		fill_pt_e( string("met0_all_slimmed_pt"), met.pt(), weight);



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

		// ---------------------------------- electrons selection
		LorentzVector elDiff(0., 0., 0., 0.);
		// std::vector<patUtils::GenericLepton>
		pat::ElectronCollection selElectrons;
		unsigned int nVetoE(0);

		for(unsigned int count_idiso_electrons = 0, n=0; n<electrons.size (); ++n)
			{
			pat::Electron& electron = electrons[n];

			bool 
				passKin(true),     passId(true),     passIso(true),
				passVetoKin(true), passVetoId(true), passVetoIso(true);

			int lid(electron.pdgId()); // should always be 11

			//apply electron corrections
			/* no lepcorrs in 13.6
			if(abs(lid)==11)
				{
				elDiff -= electron.p4();
				ElectronEnCorrector.calibrate(electron, ev.eventAuxiliary().run(), edm::StreamID::invalidStreamID()); 
				//electron = patUtils::GenericLepton(electron.el); //recreate the generic electron to be sure that the p4 is ok
				elDiff += electron.p4();
				}
			*/

			fill_pt_e( string("all_electrons_corrected_pt"), electron.pt(), weight);
			if (n < 2)
				{
				fill_pt_e( string("top2pt_electrons_corrected_pt"), electron.pt(), weight);
				}

			//no need for charge info any longer
			//lid = abs(lid);
			//TString lepStr(lid == 13 ? "mu" : "e");
			// should always be 11!
					
			// no need to mess with photon ID // //veto nearby photon (loose electrons are many times photons...)
			// no need to mess with photon ID // double minDRlg(9999.);
			// no need to mess with photon ID // for(size_t ipho=0; ipho<selPhotons.size(); ipho++)
			// no need to mess with photon ID //   minDRlg=TMath::Min(minDRlg,deltaR(leptons[n].p4(),selPhotons[ipho].p4()));
			// no need to mess with photon ID // if(minDRlg<0.1) continue;


			// ------------------------- electron IDs
			//Cut based identification
			passId     = patUtils::passId(electron, goodPV, patUtils::llvvElecId::Tight);
			passVetoId = patUtils::passId(electron, goodPV, patUtils::llvvElecId::Loose);

			// ------------------------- electron isolation

			passIso     = patUtils::passIso(electron, patUtils::llvvElecIso::Tight);
			passVetoIso = patUtils::passIso(electron, patUtils::llvvElecIso::Loose);


			// ---------------------------- kinematics
			//double leta(fabs(lid==11 ? lepton.el.superCluster()->eta() : lepton.eta()));
			double leta( electron.superCluster()->eta() );

			// ---------------------- Main lepton kin
			if(electron.pt() < 30.)                passKin = false;
			// if(leta > 2.1)                         passKin = false;
			if(leta > 2.4)                         passKin = false;
			if(leta > 1.4442 && leta < 1.5660)     passKin = false; // Crack veto

			// ---------------------- Veto lepton kin
			// if (electron.pt () < 20)            passVetoKin = false;
			if (electron.pt () < 15)            passVetoKin = false;
			// if (leta > 2.1)                     passVetoKin = false;
			if (leta > 2.5)                     passVetoKin = false;
			if (leta > 1.4442 && leta < 1.5660) passVetoKin = false; // Crack veto


			if     (passKin     && passId     && passIso)     selElectrons.push_back(electron);
			else if(passVetoKin && passVetoId && passVetoIso) nVetoE++;

			// TODO: ------------ ID/Iso scale factors:
			// https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2#Efficiencies_and_scale_factors
			// if it is MC and ID tight electron is found -- multiply weight by factor,
			// corresponding to its' tight scale factor distr
			// if loose -- loose scale factor

			/* TODO: add ID/Iso and run
			if (isMC) {
				if (passId) {
					// NOTE: scale factors are given for absolute eta
					weight *= eID_Tight_sf.GetBinContent( eID_Tight_sf.FindBin(fabs(leta), electron.pt()) );
				} else if (passVetoId) {
					weight *= eID_Loose_sf.GetBinContent( eID_Loose_sf.FindBin(fabs(leta), electron.pt()) );
				}
				if (passIso) {
					// NOTE: scale factors are given for absolute eta
					weight *= eIso_Tight_sf.GetBinContent( eIso_Tight_sf.FindBin(fabs(leta), electron.pt()) );
				} else if (passVetoIso) {
					weight *= eIso_Loose_sf.GetBinContent( eIso_Loose_sf.FindBin(fabs(leta), electron.pt()) );
				}
			}
			*/


			if (passId && passIso)
				{
				fill_pt_e( string("all_electrons_idiso_pt"), electron.pt(), weight);
				if (count_idiso_electrons < 2) fill_pt_e( string("top2pt_electrons_idiso_pt"), electron.pt(), weight);
				count_idiso_electrons += 1;
				}
			}

		// TODO: there should be no need to sort selected electrons here again -- they are in order of Pt
		std::sort(selElectrons.begin(),   selElectrons.end(),   utils::sort_CandidatesByPt);
		// std::sort(selLeptons.begin(),   selLeptons.end(),   utils::sort_CandidatesByPt);
		// std::sort(selLeptons_nocor.begin(),   selLeptons_nocor.end(),   utils::sort_CandidatesByPt);

		for(size_t n=0; n<selElectrons.size(); ++n)
			{
			fill_pt_e( string("all_electrons_pt_individual"), selElectrons[n].pt(), weight);
			if (n < 1)
				{
				fill_pt_e( string("top2pt_electrons_pt_individual"), selElectrons[n].pt(), weight);
				}
			}


		if(debug){
			cout << "processed electrons" << endl;
			}





		// ---------------------------------- muons selection
		LorentzVector muDiff(0., 0., 0., 0.);
		// std::vector<patUtils::GenericLepton> selLeptons;
		pat::MuonCollection selMuons;
		unsigned int nVetoMu(0);
		// unsigned int count_idiso_muons = 0;

		for(unsigned int count_idiso_muons = 0, n=0; n<muons.size (); ++n)
			{
			// patUtils::GenericLepton& lepton = leptons[n];
			pat::Muon& muon = muons[n];

			bool
				passKin(true),     passId(true),     passIso(true),
				passVetoKin(true), passVetoId(true), passVetoIso(true);

			int lid(muon.pdgId()); // should always be 13

			//apply muon corrections
			/* no lepcorrs in 13.6
			if(abs(lid) == 13 && muCor)
				{
				float qter;
				TLorentzVector p4(muon.px(), muon.py(), muon.pz(), muon.energy());
				// old corrections:
				// muCor->applyPtCorrection(p4, lid < 0 ? -1 : 1);
				// if(isMC) muCor->applyPtSmearing(p4, lid < 0 ? -1 : 1, false);
				// roch-cor (rochcor) corrections:
				if (isMC) muCor->momcor_mc  (p4, lid<0 ? -1 :1, 0, qter);
				else muCor->momcor_data(p4, lid<0 ? -1 :1, 0, qter);
				muDiff -= muon.p4();
				muon.setP4(LorentzVector(p4.Px(), p4.Py(), p4.Pz(), p4.E()));
				muDiff += muon.p4();
				}
			*/

			fill_pt_e( string("all_muons_corrected_pt"), muon.pt(), weight);
			if (n < 2)
				{
				fill_pt_e( string("top2pt_muons_corrected_pt"), muon.pt(), weight);
				}

			//no need for charge info any longer
			//lid = abs(lid);
			//TString lepStr(lid == 13 ? "mu" : "e");
					
			// no need to mess with photon ID // //veto nearby photon (loose electrons are many times photons...)
			// no need to mess with photon ID // double minDRlg(9999.);
			// no need to mess with photon ID // for(size_t ipho=0; ipho<selPhotons.size(); ipho++)
			// no need to mess with photon ID //   minDRlg=TMath::Min(minDRlg,deltaR(leptons[ilep].p4(),selPhotons[ipho].p4()));
			// no need to mess with photon ID // if(minDRlg<0.1) continue;

			//Cut based identification

			// ------------------------- muon IDs
			passId     = patUtils::passId(muon, goodPV, patUtils::llvvMuonId::StdTight);
			passVetoId = patUtils::passId(muon, goodPV, patUtils::llvvMuonId::StdLoose);

			// ------------------------- muon isolation
			passIso     = patUtils::passIso(muon, patUtils::llvvMuonIso::Tight);
			passVetoIso = patUtils::passIso(muon, patUtils::llvvMuonIso::Loose);

			if (passId && passIso)
				{
				fill_pt_e( string("all_muons_idiso_pt"), muon.pt(), weight);
				if (count_idiso_muons < 2) fill_pt_e( string("top2pt_muons_idiso_pt"), muon.pt(), weight);
				count_idiso_muons += 1;
				}


			// ---------------------------- kinematics
			double leta( muon.eta());

			// ---------------------- Main muon kin
			// if(muon.pt() < 30.)   passKin = false;
			// if(leta > 2.1)        passKin = false;
			if(muon.pt() < 26.)   passKin = false;
			if(leta > 2.4)        passKin = false;

			// ---------------------- Veto muon kin
			// if (muon.pt () < 20)  passVetoKin = false;
			// if (leta > 2.1)       passVetoKin = false;
			if (muon.pt () < 10)  passVetoKin = false;
			if (leta > 2.5)       passVetoKin = false;

			if     (passKin     && passId     && passIso)     selMuons.push_back(muon);
			else if(passVetoKin && passVetoId && passVetoIso) nVetoMu++;

			// TODO: ------------ ID/Iso scale factors:
			// https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffsRun2#Muon_reconstruction_identificati
			// if it is MC and ID tight muon is found -- multiply weight by factor,
			// corresponding to its' tight scale factor distr
			// if loose -- loose scale factor

			/* TODO: add ID/Iso and run
			if (isMC) {
				if (passId) {
					// NOTE: scale factors are given for absolute eta
					weight *= muID_Tight_sf.GetBinContent( muID_Tight_sf.FindBin(leta, muon.pt()) );
				} else if (passVetoId) {
					weight *= muID_Loose_sf.GetBinContent( muID_Loose_sf.FindBin(leta, muon.pt()) );
				}
				if (passIso) {
					// NOTE: scale factors are given for absolute eta
					weight *= muIso_Tight_sf.GetBinContent( muIso_Tight_sf.FindBin(leta, muon.pt()) );
				} else if (passVetoIso) {
					weight *= muIso_Loose_sf.GetBinContent( muIso_Loose_sf.FindBin(leta, muon.pt()) );
				}
			}
			*/
			}


		// TODO: there should be no need to sort selected muons here again -- they are in order of Pt
		std::sort(selMuons.begin(),   selMuons.end(),   utils::sort_CandidatesByPt);

		for(size_t n=0; n<selMuons.size(); ++n)
			{
			fill_pt_e( string("all_muons_pt_individual"), selMuons[n].pt(), weight);
			if (n < 1)
				{
				fill_pt_e( string("top2pt_muons_pt_individual"), selMuons[n].pt(), weight);
				}
			}

		if(debug){
			cout << "processed muons" << endl;
			}


		// Check for the single-electron/single-muon channels and save Pt-s if the event is in channel

		// Here we can already assign leptonic channel
		// electron or muon -- tau
		// the selection is:

		// unsigned int n_leptons = selLeptons.size();
		// Event classification. Single lepton triggers are used for offline selection of dilepton events. The "else if"s guarantee orthogonality
		bool 
			isSingleMu(false),
			isSingleE(false),
			isDoubleMu(false),
			isDoubleE(false),
			isEMu(false);
		// int multiChannel(0);


		// int slepId(0);

		// if(selLeptons.size()>0)
			// slepId=selLeptons[0].pdgId();

		// bool iso_lep = nVetoE==0 && nVetoMu==0 && selLeptons.size() == 1 && nGoodPV != 0; // 2^5
		//if(selLeptons.size()!=1 || nGoodPV==0) continue; // Veto requirement alredy applied during the event categoriziation
		// isSingleMu = (abs(slepId)==13) && muTrigger && iso_lep;
		// isSingleE  = (abs(slepId)==11) && eTrigger  && iso_lep;
		bool clean_lep_conditions = nVetoE==0 && nVetoMu==0 && nGoodPV != 0;
		// TODO: test if trigger needed here at all
		//isSingleMu = selMuons.size() == 1 && selElectrons.size() == 0 && muTrigger && clean_lep_conditions;
		//isSingleE  = selMuons.size() == 0 && selElectrons.size() == 1 && eTrigger  && clean_lep_conditions;
		// FIXME: TESTING new, trigger-less channel assignment
		isSingleMu = selMuons.size() == 1 && selElectrons.size() == 0 && clean_lep_conditions;
		isSingleE  = selMuons.size() == 0 && selElectrons.size() == 1 && clean_lep_conditions;


		if (isSingleMu) fill_pt_e( string("singlemu_muons_pt"),     selMuons[0].pt(), weight);
		if (isSingleE)  fill_pt_e( string("singleel_electrons_pt"), selElectrons[0].pt(), weight);

		if(debug){
			cout << "assigned lepton channel" << endl;
			}

		if (isSingleE)  fill_pt_e( string("met0_singleel_slimmed_pt"), met.pt(), weight);
		if (isSingleMu) fill_pt_e( string("met0_singlemu_slimmed_pt"), met.pt(), weight);

		// ------------------------------------- Propagate lepton energy scale to MET
		// Propagate now (v13)
		/* no lepton corrections propagation v13.1
		met.setP4(met.p4() - muDiff - elDiff); // TODO: note this also propagates to all MET uncertainties -- does it??
		met.setUncShift(met.px() - muDiff.px()*0.01, met.py() - muDiff.py()*0.01, met.sumEt() - muDiff.pt()*0.01, pat::MET::METUncertainty::MuonEnUp);   //assume 1% uncertainty on muon rochester
		met.setUncShift(met.px() + muDiff.px()*0.01, met.py() + muDiff.py()*0.01, met.sumEt() + muDiff.pt()*0.01, pat::MET::METUncertainty::MuonEnDown); //assume 1% uncertainty on muon rochester
		met.setUncShift(met.px() - elDiff.px()*0.01, met.py() - elDiff.py()*0.01, met.sumEt() - elDiff.pt()*0.01, pat::MET::METUncertainty::ElectronEnUp);   //assume 1% uncertainty on electron scale correction
		met.setUncShift(met.px() + elDiff.px()*0.01, met.py() + elDiff.py()*0.01, met.sumEt() + elDiff.pt()*0.01, pat::MET::METUncertainty::ElectronEnDown); //assume 1% uncertainty on electron scale correction
		*/

		fill_pt_e( string("met0_all_leptoncorr_pt"), met.pt(), weight);

		if (isSingleE)  fill_pt_e( string("met0_singleel_leptoncorr_pt"), met.pt(), weight);
		if (isSingleMu) fill_pt_e( string("met0_singlemu_leptoncorr_pt"), met.pt(), weight);

		if(debug){
			cout << "propagated lepton corrections to met" << endl;
			}

		// FIXME: this is absolutely a test procedure for leptons mismatch, delete later
		/*
		if (isSingleE && isMC)
			{
			// the ratio table then:
			double pt = selLeptons[0].pt();
			if (pt > 30 && pt < 81)
				{
				double ratios[102] = {1.417432, 1.355603, 1.326855, 1.295397, 1.274588,
					1.244550, 1.224524, 1.189812, 1.172887, 1.167444, 1.155932,
					1.126863, 1.123598, 1.101529, 1.102847, 1.102109, 1.077481,
					1.063031, 1.070628, 1.065650, 1.055504, 1.057805, 1.067119,
					1.069570, 1.073800, 1.073553, 1.090602, 1.088997, 1.089411,
					1.073635, 1.101547, 1.107392, 1.079771, 1.092838, 1.101382,
					1.110520, 1.104568, 1.106222, 1.091441, 1.105385, 1.087636,
					1.090410, 1.123266, 1.116482, 1.159609, 1.112731, 1.140635,
					1.193607, 1.153018, 1.111796, 1.169275, 1.151467, 1.159238,
					1.203220, 1.126885, 1.160726, 1.140825, 1.133395, 1.171352,
					1.115036, 1.197012, 1.210529, 1.143451, 1.229832, 1.207743,
					1.202648, 1.239188, 1.228262, 1.170988, 1.289559, 1.133405,
					1.201752, 1.173881, 1.273580, 1.159010, 1.163254, 1.181294, 1.128886,
					1.260505, 1.217777, 1.113412, 1.177267, 1.230890, 1.250940,
					1.261494, 1.234183, 1.187044, 1.205113, 1.248632, 1.280663, 1.216726,
					1.225326, 1.291109, 1.256174, 1.199306, 1.266343, 1.235425,
					1.281472, 1.219297, 1.318849, 1.164462, 1.416536};
				weight *= ratios[int((pt - 30)*2)];
				}
			}

		if (isSingleMu && isMC)
			{
			double pt = selLeptons[0].pt();
			if (pt > 30 && pt < 45.5)
				{
				double ratios[31] = {1.163884, 1.150239, 1.132176, 1.113742, 1.096949,
					1.078991, 1.082796, 1.077936, 1.065941, 1.065551, 1.065027,
					1.069034, 1.058150, 1.053311, 1.052812, 1.044186, 1.036055,
					1.041379, 1.039200, 1.035185, 1.034690, 1.032413, 1.025401,
					1.022467, 1.026297, 1.023207, 1.026607, 1.023327, 1.003550,
					1.017521, 1.026832};
				weight *= ratios[int((pt - 30)*2)];
				}
			}
		*/

		// Finally, merge leptons for cross-cleaning with taus and jets:

		std::vector<patUtils::GenericLepton> selLeptons;
		for(size_t l=0; l<selElectrons.size(); ++l) selLeptons.push_back(patUtils::GenericLepton (selElectrons[l] ));
		for(size_t l=0; l<selMuons.size(); ++l)     selLeptons.push_back(patUtils::GenericLepton (selMuons[l]     ));
		std::sort(selLeptons.begin(), selLeptons.end(), utils::sort_CandidatesByPt);


		if(debug){
			cout << "merged selected leptons" << endl;
			}

		// TODO: more conditions for double-lepton channel? No veto leptons etc?
		// in progress...
		if (selLeptons.size()==2 && clean_lep_conditions)
			{
			// so this should be a double-lepton channel
			// TODO: remove tau-selection below?
			int dilep_ids = selLeptons[0].pdgId() * selLeptons[1].pdgId();

			if (fabs(dilep_ids) == 121 )
				{
				isDoubleE = true;
				fill_pt_e( string("leptons_doublee_2leptons_pt"), selLeptons[0].pt(), weight);
				fill_pt_e( string("leptons_doublee_2leptons_pt"), selLeptons[1].pt(), weight);
				}
			else if (fabs(dilep_ids) == 169 )
				{
				isDoubleMu = true;
				fill_pt_e( string("leptons_doublemu_2leptons_pt"), selLeptons[0].pt(), weight);
				fill_pt_e( string("leptons_doublemu_2leptons_pt"), selLeptons[1].pt(), weight);
				}
			else
				{
				isEMu = true;
				fill_pt_e( string("leptons_emu_2leptons_pt"), selLeptons[0].pt(), weight);
				fill_pt_e( string("leptons_emu_2leptons_pt"), selLeptons[1].pt(), weight);
				}
			}

		if (isSingleE || isSingleMu || isDoubleMu || isDoubleE || isEMu ){
			fill_pu( string("pileup_inachannel_rawweight_pergoodpv"), nGoodPV, rawWeight);
			fill_pu( string("pileup_inachannel_weight_pergoodpv"), nGoodPV, weight);
			fill_pu( string("pileup_inachannel_weight_up_pergoodpv"), nGoodPV, weight_up);
			fill_pu( string("pileup_inachannel_weight_down_pergoodpv"), nGoodPV, weight_down);

			fill_pu( string("pileup_inachannel_rawweight_perrawvtxsize"), vtx.size(), rawWeight);
			fill_pu( string("pileup_inachannel_weight_perrawvtxsize"), vtx.size(), weight_pu_test);

			fill_pu( string("pileup_inachannel_rawweight_pernuminters"), num_inters, rawWeight);
			fill_pu( string("pileup_inachannel_weight_pernuminters"), num_inters, weight);
			fill_pu( string("pileup_inachannel_weight_up_pernuminters"), num_inters, weight_up);
			fill_pu( string("pileup_inachannel_weight_down_pernuminters"), num_inters, weight_down);
			}

		/* old lepton selection, left for reference
		// ---------------------------------- leptons selection
		LorentzVector muDiff(0., 0., 0., 0.), elDiff(0., 0., 0., 0.);
		std::vector<patUtils::GenericLepton> selLeptons;
		std::vector<patUtils::GenericLepton> selLeptons_nocor;
		unsigned int nVetoE(0), nVetoMu(0);
		unsigned int nVetoE_nocor(0), nVetoMu_nocor(0);
		for(size_t ilep=0; ilep<leptons.size (); ++ilep)
			{
			patUtils::GenericLepton& lepton = leptons[ilep];
			patUtils::GenericLepton  nocor_lepton = leptons[ilep];

			bool 
				passKin(true),     passId(true),     passIso(true),
				passVetoKin(true), passVetoId(true), passVetoIso(true);

			int lid(lepton.pdgId());

			//apply muon corrections
			if(abs(lid) == 13 && muCor)
				{
				float qter;
				TLorentzVector p4(lepton.px(), lepton.py(), lepton.pz(), lepton.energy());
				// old corrections:
				// muCor->applyPtCorrection(p4, lid < 0 ? -1 : 1);
				// if(isMC) muCor->applyPtSmearing(p4, lid < 0 ? -1 : 1, false);
				// roch-cor (rochcor) corrections:
				if (isMC) muCor->momcor_mc  (p4, lid<0 ? -1 :1, 0, qter);
				else muCor->momcor_data(p4, lid<0 ? -1 :1, 0, qter);
				muDiff -= lepton.p4();
				lepton.setP4(LorentzVector(p4.Px(), p4.Py(), p4.Pz(), p4.E()));
				muDiff += lepton.p4();
				}

			//apply electron corrections
			if(abs(lid)==11)
				{
				elDiff -= lepton.p4();
				ElectronEnCorrector.calibrate(lepton.el, ev.eventAuxiliary().run(), edm::StreamID::invalidStreamID()); 
				lepton = patUtils::GenericLepton(lepton.el); //recreate the generic lepton to be sure that the p4 is ok
				elDiff += lepton.p4();
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
			// passId     = lid == 11 ? patUtils::passId(electronVidMainId, myEvent, lepton.el) : patUtils::passId(lepton.mu, goodPV, patUtils::llvvMuonId::StdTight);
			// passVetoId = lid == 11 ? patUtils::passId(electronVidVetoId, myEvent, lepton.el) : patUtils::passId(lepton.mu, goodPV, patUtils::llvvMuonId::StdLoose);
			// apparently the previous version of the passID callse is actually incompatible with the definition of passID
			// don't know how it compiled at all....
			passId     = lid == 11 ? patUtils::passId(lepton.el, goodPV, patUtils::llvvElecId::Tight) : patUtils::passId(lepton.mu, goodPV, patUtils::llvvMuonId::StdTight);
			passVetoId = lid == 11 ? patUtils::passId(lepton.el, goodPV, patUtils::llvvElecId::Loose) : patUtils::passId(lepton.mu, goodPV, patUtils::llvvMuonId::StdLoose);

			// ------------------------- lepton isolation
			// passIso     = lid == 11 ? true : patUtils::passIso(lepton.mu, patUtils::llvvMuonIso::Tight); // Electron iso is included within the ID
			// passVetoIso = lid == 11 ? true : patUtils::passIso(lepton.mu, patUtils::llvvMuonIso::Loose); // Electron iso is included within the ID
			passIso     = lid == 11 ? patUtils::passIso(lepton.el, patUtils::llvvElecIso::Tight) : patUtils::passIso(lepton.mu, patUtils::llvvMuonIso::Tight);
			passVetoIso = lid == 11 ? patUtils::passIso(lepton.el, patUtils::llvvElecIso::Loose) : patUtils::passIso(lepton.mu, patUtils::llvvMuonIso::Loose);

			if     (passKin     && passId     && passIso)     selLeptons.push_back(lepton);
			else if(passVetoKin && passVetoId && passVetoIso) lid==11 ? nVetoE++ : nVetoMu++;



			// -----------  also, store the Pt distr for passing but not corrected lepton: {TEST}

			passKin     = true;
			passId      = true;
			passIso     = true;
			passVetoKin = true;
			passVetoId  = true;
			passVetoIso = true;

			leta = fabs(lid==11 ? nocor_lepton.el.superCluster()->eta() : nocor_lepton.eta());
			
			// ---------------------- Main lepton kin
			if(nocor_lepton.pt() < 30.)                       passKin = false;
			if(leta > 2.1)                                    passKin = false;
			if(lid == 11 && (leta > 1.4442 && leta < 1.5660)) passKin = false; // Crack veto

			// ---------------------- Veto lepton kin
			if (nocor_lepton.pt () < 20)                       passVetoKin = false;
			if (leta > 2.1)                                    passVetoKin = false;
			if (lid == 11 && (leta > 1.4442 && leta < 1.5660)) passVetoKin = false; // Crack veto

			// ------------------------- lepton IDs
			passId     = lid == 11 ? patUtils::passId(nocor_lepton.el, goodPV, patUtils::llvvElecId::Tight) : patUtils::passId(nocor_lepton.mu, goodPV, patUtils::llvvMuonId::StdTight);
			passVetoId = lid == 11 ? patUtils::passId(nocor_lepton.el, goodPV, patUtils::llvvElecId::Loose) : patUtils::passId(nocor_lepton.mu, goodPV, patUtils::llvvMuonId::StdLoose);

			// ------------------------- lepton isolation
			passIso     = lid == 11 ? patUtils::passIso(nocor_lepton.el, patUtils::llvvElecIso::Tight) : patUtils::passIso(nocor_lepton.mu, patUtils::llvvMuonIso::Tight);
			passVetoIso = lid == 11 ? patUtils::passIso(nocor_lepton.el, patUtils::llvvElecIso::Loose) : patUtils::passIso(nocor_lepton.mu, patUtils::llvvMuonIso::Loose);

			if     (passKin     && passId     && passIso)     selLeptons_nocor.push_back(nocor_lepton);
			else if(passVetoKin && passVetoId && passVetoIso) lid==11 ? nVetoE_nocor++ : nVetoMu_nocor++;
			}

		std::sort(selLeptons.begin(),   selLeptons.end(),   utils::sort_CandidatesByPt);
		std::sort(selLeptons_nocor.begin(),   selLeptons_nocor.end(),   utils::sort_CandidatesByPt);
		//LorentzVector recoMET = met;// FIXME REACTIVATE IT - muDiff;

		*/


		// ------------------------------------------ select the individual taus
		pat::TauCollection selTaus;
		int ntaus (0);
		for (unsigned int count_ided_taus = 0, n = 0; n < taus.size(); ++n)
			{
			pat::Tau& tau = taus[n];

			// ---------- IDs
					
			// if(!tau.isPFTau()) continue; // Only PFTaus // It should be false for slimmedTaus
			// if(tau.emFraction() >=2.) continue;
					
			// Discriminators from https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV
			// "The tau passes the discriminator if pat::Tau::tauID("name") returns a value of 0.5 or greater"
			//if(tau.tauID("decayModeFindingNewDMs")<0.5) continue; // High pt tau. Otherwise, OldDMs
			if(tau.tauID("decayModeFinding")<0.5) continue; // High pt tau. Otherwise, OldDMs (or no <DMs> -- they are synonyms)
			// Anyways, the collection of taus from miniAOD should be already afer decayModeFinding cut (the tag - Old or New - is unspecified in the twiki, though).
			// Consequently, there might be a small bias due to events that are cut by the OldDM and would not be cut by the NewDM
			//if (tau.tauID ("byMediumCombinedIsolationDeltaBetaCorr3Hits")<0.5) continue; // See whether to us the new byMediumPileupWeightedIsolation3Hits that is available only for dynamic strip reconstruction (default in CMSSW_7_4_14)
			if (tau.tauID ("byTightCombinedIsolationDeltaBetaCorr3Hits")<0.5) continue;
			if (tau.tauID ("againstMuonTight3")                          <0.5) continue; // Medium working point not usable. Available values: Loose, Tight
			//if (tau.tauID ("againstElectronMediumMVA5")                  <0.5) continue; // Tight working point not usable. Avaiable values: VLoose, Loose, Medium
			// if (tau.tauID ("againstElectronMediumMVA6")                  <0.5) continue;
			if (tau.tauID ("againstElectronTightMVA6")                  <0.5) continue;
			
			fill_pt_e( string("all_taus_4discrs_pt"), tau.pt(), weight);

			// Pixel hits cut (will be available out of the box in new MINIAOD production)
			{
			int nChHadPixelHits = 0;
			reco::CandidatePtrVector chCands = tau.signalChargedHadrCands();
			for(reco::CandidatePtrVector::const_iterator iter = chCands.begin(); iter != chCands.end(); iter++)
				{
				pat::PackedCandidate const* packedCand = dynamic_cast<pat::PackedCandidate const*>(iter->get());
				int pixelHits = packedCand->numberOfPixelHits();
				if(pixelHits > nChHadPixelHits) nChHadPixelHits = pixelHits;
				}
			if(nChHadPixelHits==0) continue;
			}

			fill_pt_e( string("all_taus_ided_pt"), tau.pt(), weight);
			if (count_ided_taus<1)
				{
				fill_pt_e( string("top1pt_taus_ided_pt"), tau.pt(), weight);
				count_ided_taus += 1;
				}

			// --------- Kinematics
			if (tau.pt() < 20. || fabs (tau.eta()) > 2.3) continue;

			selTaus.push_back(tau);
			ntaus++;
			}
		std::sort (selTaus.begin(), selTaus.end(), utils::sort_CandidatesByPt);

		// Control values for processed individual taus:
		for(size_t n=0; n<selTaus.size(); ++n)
			{
			fill_pt_e( string("all_taus_individual_pt"), selTaus[n].pt(), weight);
			if (n == 0)
				{
				fill_pt_e( string("top1pt_taus_individual_pt"), selTaus[n].pt(), weight);
				if (isSingleE)  fill_pt_e( string("top1pt_taus_individual_singleel_pt"), selTaus[n].pt(), weight);
				if (isSingleMu) fill_pt_e( string("top1pt_taus_individual_singlemu_pt"), selTaus[n].pt(), weight);
				}
			}

		if(debug){
			cout << "selected taus [individual]" << endl;
			}


		// ------------------------------------------ select the taus cleaned from leptons

		pat::TauCollection selTausNoLep;
		int closest_totaunolep_particle_id = 0; // wonder what is 0 particle
		for (size_t itau = 0; itau < selTaus.size(); ++itau)
			{
			pat::Tau& tau = selTaus[itau];

			// cross-cleaning taus with leptons
			bool overlapWithLepton(false);
			for(int l=0; l<(int)selLeptons.size();++l)
				{
				if (reco::deltaR(tau, selLeptons[l])<0.4)
					{ overlapWithLepton=true; break; }
				}
			if (overlapWithLepton) continue;

			selTausNoLep.push_back(tau);
			// so these are the final taus we use in the selection

			// for the fake-rate counts (in MC)
			// let's save how many taus we find:
			increment(string("number_of_tausnolep_found"), 1);

			// for MC find the generated candidate closest to the tau
			// to get the fakeness
			if (isMC)
				{
				double min_deltaR = 99999.9;
				for(size_t i = 0; i < gen.size(); ++ i)
					{
					const reco::GenParticle & p = gen[i];
					double deltaR_to_p = reco::deltaR(tau, p);
					if (deltaR_to_p < min_deltaR)
						{
						min_deltaR = deltaR_to_p;
						closest_totaunolep_particle_id = p.pdgId();
						}
					//int id = p.pdgId();
					//int st = p.status();
					//int n = p.numberOfDaughters();
					//cout << i << ": " << id << " " << st;
					//if (p.numberOfMothers() != 0) cout << " <- " ;
					//for (int j = 0 ; j < p.numberOfMothers(); ++j)
					//	{
					//	const reco::Candidate * mom = p.mother(j);
					//	cout << " " << mom->pdgId() << " " << mom->status() << ";";
					//	}
					//cout << "\n";
					//if (n>0)
					//	{
					//	cout << "\t|-> " ;
					//	for (int j = 0; j < n; ++j)
					//		{
					//		const reco::Candidate * d = p.daughter( j );
					//		cout << d->pdgId() << " " << d->status() << "; " ;
					//		}
					//	cout << "\n";
					//	}
					}
				fill_particle_ids(string("nearest_particle_around_taunolep"), closest_totaunolep_particle_id, weight);
				// electron-tau fake-rate scale factor
				if (fabs(closest_totaunolep_particle_id)==11)
					{
					increment(string("number_of_tausnolep_from_electron_found"), 1);
					// apply the data/MC fake rate scale factor
					// (from https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV#Electron_to_tau_fake_rate)
					// for tight working points it is 1.80 +- 0.23 in barrel (eta < 1.460)
					// and 1.3 +- 0.42 in endcaps (eta > 1.558)
					/* repeating 13.4, no electron-to-tau sf
					if (fabs(tau.eta()) < 1.460)
						{
						weight *= 1.8 + r3->Gaus(0, 0.2);  // gaussian +- 0.2
						}
					else if (fabs(tau.eta()) > 1.558)
						{
						weight *= 1.3 + r3->Gaus(0, 0.42); // gaussian +- 0.42
						}
					*/
					// TODO: and then I need to renormalize the MC integral??
					}

				if (fabs(closest_totaunolep_particle_id)==13)
					{
					increment(string("number_of_tausnolep_from_muon_found"), 1);
					// TODO: add muon-tau fake rate?
					}
				}
			}

		increment(string("weightflow_weight_after_tausnolep_fakerates_sf"), weight);

		// Control values for processed taus cleaned of leptons:
		for(size_t n=0; n<selTausNoLep.size(); ++n)
			{
			fill_pt_e( string("all_taus_leptoncleaned_pt"), selTausNoLep[n].pt(), weight);
			if (n == 0)
				{
				fill_pt_e( string("top1pt_taus_leptoncleaned_pt"), selTausNoLep[n].pt(), weight);
				if (isSingleE)  fill_pt_e( string("top1pt_taus_lepcleaned_singleel_pt"), selTausNoLep[n].pt(), weight);
				if (isSingleMu) fill_pt_e( string("top1pt_taus_lepcleaned_singlemu_pt"), selTausNoLep[n].pt(), weight);
				}
			}



		if(debug){
			cout << "processed taus" << endl;
			}

		//
		// ----------------------------------------------- JET/MET ANALYSIS
		//
		if(debug) cout << "Now update Jet Energy Corrections" << endl;
		//add scale/resolution uncertainties and propagate to the MET
		//utils::cmssw::updateJEC(jets, jesCor, totalJESUnc, rho, nGoodPV, isMC);

		// up to here jets were not processed in any way
		// now goes the procedure of corrections to jets and then METs

		// ----------------------- The updateJEC procedure from src/MacroUtils:
		//void updateJEC(pat::JetCollection& jets, FactorizedJetCorrector *jesCor, JetCorrectionUncertainty *totalJESUnc, float rho, int nvtx,bool isMC){

		// inline control functions usage:
		//   fill_pt_e( "control_point_name", value, weight)
		//   fill_eta( "control_point_name", value, weight )   <-- different TH1F range and binning
		//   increment( "control_point_name", weight )


		if(debug) cout << "jet eta pt e, e x y z" << endl;

		LorentzVector jet_corr(0., 0., 0., 0.);
		for(size_t ijet=0; ijet<jets.size(); ijet++)
			{
			// TODO: so does this mean "in place"?
			pat::Jet& jet = jets[ijet];
			LorentzVector jet_initial_momentum = jet.p4();

			if(debug) cout << jet.eta() << " " << jet.pt() << " " << jet.energy() << endl;
			//correct JES
			LorentzVector rawJet = jet.correctedP4("Uncorrected");

			if(debug) cout << rawJet.eta() << " " << rawJet.pt() << " " << rawJet.energy() << endl;

			// here is the correct1 jet correction point

			fill_pt_e( string("all_jets_pt_raw"), rawJet.pt(), weight);
			if (ijet < 2)
				{
				fill_pt_e( string("top2pt_jets_pt_raw"), rawJet.pt(), weight);
				}

			//double toRawSF=jet.correctedJet("Uncorrected").pt()/jet.pt();
			//LorentzVector rawJet(jet*toRawSF);
			// 13.8_2 trying all jet corrections + new (if it is new) jetID procedure
			jesCor->setJetEta(rawJet.eta());
			jesCor->setJetPt(rawJet.pt());
			jesCor->setJetA(jet.jetArea());
			jesCor->setRho(rho);
			jesCor->setNPV(nGoodPV);
			jet.setP4(rawJet*jesCor->getCorrection());

			// here is the correct2 jet correction point

			fill_pt_e( string("all_jets_pt_corrected2"), jet.pt(), weight);
			if (ijet < 2)
				{
				fill_pt_e( string("top2pt_jets_pt_corrected2"), jet.pt(), weight);
				}


			if(debug) cout << jet.eta() << " " << jet.pt() << " " << jet.energy() << endl;

			//smear JER
			//double newJERSF(1.0);
			// 13.8_2 trying all jet corrections + new (if it is new) jetID procedure
			if(isMC)
				{
				const reco::GenJet* genJet=jet.genJet();
				if(genJet)
					{
					double genjetpt( genJet ? genJet->pt(): 0.);                    
					std::vector<double> smearJER=utils::cmssw::smearJER(jet.pt(),jet.eta(),genjetpt);
					jet.setP4(jet.p4()*smearJER[0]);
					
					//printf("jet pt=%f gen pt = %f smearing %f %f %f\n", jet.pt(), genjetpt, smearJER[0], smearJER[1], smearJER[2]);
					// //set the JER up/down alternatives
					jet.addUserFloat("jerup", smearJER[1]);  //kept for backward compatibility
					jet.addUserFloat("jerdown", smearJER[2] ); //kept for backward compatibility
					jet.addUserFloat("_res_jup", smearJER[1]);
					jet.addUserFloat("_res_jdown", smearJER[2] );
					}
				else{
					jet.addUserFloat("jerup", 1.0); //kept for backward compatibility
					jet.addUserFloat("jerdown", 1.0);  //kept for backward compatibility
					jet.addUserFloat("_res_jup", 1.0);
					jet.addUserFloat("_res_jdown", 1.0 );
					}
				}

			// here is the correct3 jet correction point

			if(isMC)
				{
				////set the JES up/down pT alternatives
				std::vector<float> ptUnc=utils::cmssw::smearJES(jet.pt(),jet.eta(), totalJESUnc);
				jet.addUserFloat("jesup",    ptUnc[0] );  //kept for backward compatibility
				jet.addUserFloat("jesdown",  ptUnc[1] );  //kept for backward compatibility
				jet.addUserFloat("_scale_jup",    ptUnc[0] );
				jet.addUserFloat("_scale_jdown",  ptUnc[1] );
				}

			// FIXME: this is not to be re-set. Check that this is a desired non-feature.
			// i.e. check that the uncorrectedJet remains the same even when the corrected momentum is changed by this routine.
			//to get the raw jet again
			//jets[ijet].setVal("torawsf",1./(newJECSF*newJERSF));

			// Add the jet momentum correction:
			// jet_cor propagation is on in 13.4
			jet_corr += jet.p4() - jet_initial_momentum;

			if(debug)
				{
				cout << jet.eta() << " " << jet.pt() << " " << jet.energy() << endl;
				cout << "-" << jet_initial_momentum.eta() << " " << jet_initial_momentum.pt() << " " << jet_initial_momentum.energy() << endl;
				}
			}


		std::sort (jets.begin(),  jets.end(),  utils::sort_CandidatesByPt);
		// ----------------------------------- here is the correctF jet correction point
		// Propagate jet_corr to MET:
		met.setP4(met.p4() - jet_corr);
		// TODO: uncertainties?
		// for leptons they are done as:
		//met.setUncShift(met.px() - muDiff.px()*0.01, met.py() - muDiff.py()*0.01, met.sumEt() - muDiff.pt()*0.01, pat::MET::METUncertainty::MuonEnUp);   //assume 1% uncertainty on muon rochester
		//met.setUncShift(met.px() + muDiff.px()*0.01, met.py() + muDiff.py()*0.01, met.sumEt() + muDiff.pt()*0.01, pat::MET::METUncertainty::MuonEnDown); //assume 1% uncertainty on muon rochester
		//met.setUncShift(met.px() - elDiff.px()*0.01, met.py() - elDiff.py()*0.01, met.sumEt() - elDiff.pt()*0.01, pat::MET::METUncertainty::ElectronEnUp);   //assume 1% uncertainty on electron scale correction
		//met.setUncShift(met.px() + elDiff.px()*0.01, met.py() + elDiff.py()*0.01, met.sumEt() + elDiff.pt()*0.01, pat::MET::METUncertainty::ElectronEnDown); //assume 1% uncertainty on electron scale correction



		// Control values for corrected jets:


		for(size_t n=0; n<jets.size(); ++n)
			{
			if(debug) cout << jets[n].eta() << " " << jets[n].pt() << " " << jets[n].energy() << endl;
			fill_pt_e( string("all_jets_pt_correctedF"), jets[n].pt(), weight);
			if (n < 2)
				{
				fill_pt_e( string("top2pt_jets_pt_correctedF"), jets[n].pt(), weight);
				}
			}



		// FIXME: So are these MET corrections?
		if(debug) cout << "Update also MET" << endl;
		// LorentzVector n_met = met.p4();
		// std::vector<LorentzVector> newMet = utils::cmssw::getMETvariations(n_met/*recoMet*/,jets,selLeptons,isMC);
		// FIXME: Must choose a lepton collection. Perhaps loose leptons?
		// n_met = newMet[utils::cmssw::METvariations::NOMINAL];

		//fill_pt_e( string("met0_all_leptoncorr_jetcorr_pt"), n_met.pt(), weight);
		//if (isSingleMu) fill_pt_e( string("met0_all_leptoncorr_jetcorr_singlemu_pt"), n_met.pt(), weight);
		//if (isSingleE)  fill_pt_e( string("met0_all_leptoncorr_jetcorr_singleel_pt"), n_met.pt(), weight);
		fill_pt_e( string("met0_all_leptoncorr_jetcorr_pt"), met.pt(), weight);
		if (isSingleMu) fill_pt_e( string("met0_all_leptoncorr_jetcorr_singlemu_pt"), met.pt(), weight);
		if (isSingleE)  fill_pt_e( string("met0_all_leptoncorr_jetcorr_singleel_pt"), met.pt(), weight);

		if(debug) cout << "Jet Energy Corrections updated" << endl;


		// TODO: should MET corrections be done here?
		// METs with corrections
		// not used now
		//
		//double met_pt_values[7];
		//met_pt_values[0] = n_met.pt();
		//met_pt_values[1] = mets[0].shiftedP4(pat::MET::METUncertainty::JetEnUp).pt();
		//met_pt_values[2] = mets[0].shiftedP4(pat::MET::METUncertainty::JetEnDown).pt();
		//met_pt_values[3] = mets[0].shiftedP4(pat::MET::METUncertainty::JetResUp).pt();
		//met_pt_values[4] = mets[0].shiftedP4(pat::MET::METUncertainty::JetResDown).pt();
		//met_pt_values[5] = mets[0].shiftedP4(pat::MET::METUncertainty::UnclusteredEnUp).pt();
		//met_pt_values[6] = mets[0].shiftedP4(pat::MET::METUncertainty::UnclusteredEnDown).pt();

		// ------------------------------- Select the jets, based on their individual parameters
		// I need different collections because of tau cleaning, but this is needed only for the single lepton channels, so the tau cleaning is performed later.
		pat::JetCollection selJets;
		// TODO: do all jet selection right here
		// now selBJets are not used anywhere
		// selJets pass cross-cleaning with taus later
		// and b-tagging again
		double mindphijmet (9999.);
		for (unsigned int count_ided_jets = 0, ijet = 0; ijet < jets.size(); ++ijet)
			{
			pat::Jet& jet = jets[ijet];

			// TODO: what do we do here exactly?
			// a loose selection on jets, and then tighten it later?
			// if (jet.pt() < 15 || fabs (jet.eta()) > 3.0) continue;
			// Was 4.7 in eta. Tightened for computing time. 3.0 ensures that we don't cut associations with leptons (0.4 from 2.4)

			//mc truth for this jet
			//const reco::GenJet * genJet = jet.genJet();
			//TString jetType (genJet && genJet->pt() > 0 ? "truejetsid" : "pujetsid");
			// TODO: this mctruth for jets it is never used in the code

			//jet id
			//bool passPFloose = passPFJetID("Loose", jet); 
			bool passPFloose = passPFJetID("Tight", jet); 
			//if (label=="Tight")
			// FIXME: check when pileup ID will come out

			if (passPFloose)
				{
				fill_pt_e( string("all_jets_ided_pt"), jet.pt(), weight);
				if (count_ided_jets < 2)
					{
					fill_pt_e( string("top2pt_jets_ided_pt"), jet.pt(), weight);
					count_ided_jets += 1;
					}
				}

			// and now the tighter final selection
			double eta = jet.eta();
			double pt  = jet.pt();
			bool passKino = pt > 30. && fabs(eta) < 2.5;

			// corrections:
			// TODO: are they MC-only?
			/* no smeared  jet pt values for now
			std::vector<double> pt_values;
			if (isMC)
				pt_values = utils::cmssw::smearJES(pt, eta, totalJESUnc);
			else
				{
				pt_values.push_back(pt);
				pt_values.push_back(pt);
				}
			*/
			// vary JesUp   is pt_values[0]
			// vary JesDown is pt_values[1]
			// if (!passPFloose || jet.pt() <30. || fabs(jet.eta()) > 2.5) continue;
			// if (passPFloose && (pt > 30. || pt_values[0] > 30. || pt_values[1] > 30.) && fabs(eta) < 2.5)


			if (passPFloose && passKino)
				{
				selJets.push_back(jet);

				// double dphijmet = fabs (deltaPhi (n_met.phi(), jet.phi()));
				double dphijmet = fabs (deltaPhi (met.phi(), jet.phi()));
				if (dphijmet < mindphijmet) mindphijmet = dphijmet;
				// FIXME: mindphijmet is not used anywhere now
				}
			}

		std::sort (selJets.begin(),  selJets.end(),  utils::sort_CandidatesByPt);

		// Control values for processed individual jets:
		for(size_t n=0; n<selJets.size(); ++n)
			{

			fill_pt_e( string("all_jets_individual_pt"), selJets[n].pt(), weight);
			if (n < 2)
				{
				fill_pt_e( string("top2pt_jets_individual_pt"), selJets[n].pt(), weight);
				}
			}




		// ---------------------------- Clean jet collections from selected leptons
		// TODO: add gamma-cleaning as well?

		pat::JetCollection selJetsNoLep;
		for (size_t ijet = 0; ijet < selJets.size(); ++ijet)
			{
			pat::Jet& jet = selJets[ijet];

			double minDRlj (9999.);

			for (size_t ilep = 0; ilep < selLeptons.size(); ilep++)
				minDRlj = TMath::Min(minDRlj, reco::deltaR (jet, selLeptons[ilep]));

			if (minDRlj < 0.4) continue;

			selJetsNoLep.push_back(jet);
			}


		// Control values for processed jets cleaned of leptons:
		for(size_t n=0; n<selJetsNoLep.size(); ++n)
			{

			fill_pt_e( string("all_jets_leptoncleaned_pt"), selJetsNoLep[n].pt(), weight);
			if (n < 2)
				{
				fill_pt_e( string("top2pt_jets_leptoncleaned_pt"), selJetsNoLep[n].pt(), weight);
				}
			}



		// ---------------------------- Clean jet collection from selected taus
		pat::JetCollection selJetsNoLepNoTau;

		for (size_t ijet = 0; ijet < selJetsNoLep.size(); ++ijet)
			{
			pat::Jet& jet = selJetsNoLep[ijet];

			double minDRtj(9999.);

			for(size_t itau=0; itau < selTausNoLep.size(); ++itau)
				minDRtj = TMath::Min(minDRtj, reco::deltaR(jet, selTausNoLep[itau]));

			if (minDRtj < 0.4) continue;

			selJetsNoLepNoTau.push_back(jet);
			}

		// Control values for processed jets cleaned of leptons and taus:
		for(size_t n=0; n<selJetsNoLepNoTau.size(); ++n)
			{

			fill_pt_e( string("all_jets_taucleaned_pt"), selJetsNoLepNoTau[n].pt(), weight);
			if (n < 2)
				{
				fill_pt_e( string("top2pt_jets_taucleaned_pt"), selJetsNoLepNoTau[n].pt(), weight);
				}

			if (isSingleMu)
				{
				fill_pt_e( string("all_jets_taucleaned_singlemu_pt"), selJetsNoLepNoTau[n].pt(), weight);
				if (n < 2)
					{
					fill_pt_e( string("top2pt_jets_taucleaned_singlemu_pt"), selJetsNoLepNoTau[n].pt(), weight);
					}
				}

			if (isSingleE)
				{
				fill_pt_e( string("all_jets_taucleaned_singleel_pt"), selJetsNoLepNoTau[n].pt(), weight);
				if (n < 2)
					{
					fill_pt_e( string("top2pt_jets_taucleaned_singleel_pt"), selJetsNoLepNoTau[n].pt(), weight);
					}
				}

			if (isEMu)
				{
				fill_pt_e( string("jets_all_taucleaned_emu_pt"), selJetsNoLepNoTau[n].pt(), weight);
				if (n < 2)
					{
					fill_pt_e( string("jets_top2pt_taucleaned_emu_pt"), selJetsNoLepNoTau[n].pt(), weight);
					}
				}
			}


		if(debug){
			cout << "processed jets" << endl;
			}


		// --------------------------- B-tagged jets
		pat::JetCollection selBJets;

		for (size_t ijet = 0; ijet < selJetsNoLepNoTau.size(); ++ijet)
			{
			pat::Jet& jet = selJetsNoLepNoTau[ijet];

			bool hasCSVtag(jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > btagMedium); // old working point
			// bool hasCSVtag(jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.8); // new working point -- according to Mara's analysis
			bool hasCSVtag_BTagUp(false), hasCSVtag_BTagDown(false);

			//update according to the SF measured by BTV
			// TODO: new fency procedure with CSV files
			if(isMC){
				// int flavId=jet.partonFlavour();
				int flavId=jet.hadronFlavour();
				double eta=jet.eta();
				if (abs(flavId)==5) {
					btsfutil.modifyBTagsWithSF(hasCSVtag, btagCal.eval(BTagEntry::FLAV_B   , eta, jet.pt()), beff);
					// btsfutil.modifyBTagsWithSF(hasCSVtagUp  , btagCalUp .eval(BTagEntry::FLAV_B   , eta, jet.pt()), beff);
					// btsfutil.modifyBTagsWithSF(hasCSVtagDown, btagCalDn .eval(BTagEntry::FLAV_B   , eta, jet.pt()), beff);
				} else if(abs(flavId)==4) {
					btsfutil.modifyBTagsWithSF(hasCSVtag, btagCal.eval(BTagEntry::FLAV_C   , eta, jet.pt()), beff);
					// btsfutil.modifyBTagsWithSF(hasCSVtagUp  , btagCalUp .eval(BTagEntry::FLAV_C   , eta, jet.pt()), beff);
					// btsfutil.modifyBTagsWithSF(hasCSVtagDown, btagCalDn .eval(BTagEntry::FLAV_C   , eta, jet.pt()), beff);
				} else {
					btsfutil.modifyBTagsWithSF(hasCSVtag, btagCalL.eval(BTagEntry::FLAV_UDSG, eta, jet.pt()), leff);
					// btsfutil.modifyBTagsWithSF(hasCSVtagUp  , btagCalLUp.eval(BTagEntry::FLAV_UDSG, eta, jet.pt()), leff);
					// btsfutil.modifyBTagsWithSF(hasCSVtagDown, btagCalLDn.eval(BTagEntry::FLAV_UDSG, eta, jet.pt()), leff);
				}
			}

			/*
			//update according to the SF measured by BTV
			if (isMC)
				{
				// test:
				//int flavId = jet.partonFlavour();
				//if      (abs(flavId)==5) btsfutil.modifyBTagsWithSF(hasCSVtag, sfb,   beff);
				//else if (abs(flavId)==4) btsfutil.modifyBTagsWithSF(hasCSVtag, sfb/5, beff);
				//else                     btsfutil.modifyBTagsWithSF(hasCSVtag, sfl,   leff);

				// TODO: also Pietro now has a more complex modifyBTagsWithSF:
				//      btsfutil.modifyBTagsWithSF(hasCSVtagDown, btagCalDn.eval(BTagEntry::FLAV_B   , eta, jet.pt()), beff);
				//  --- etc

				/* TODO: for later
				if      (abs(flavId)==5) btsfutil.modifyBTagsWithSF(hasCSVtag_BTagUp, sfb + sfbunc,   beff);
				else if (abs(flavId)==4) btsfutil.modifyBTagsWithSF(hasCSVtag_BTagUp, sfb/5 + 2*sfbunc, beff);
				else                     btsfutil.modifyBTagsWithSF(hasCSVtag_BTagUp, sfl + sfbunc,   leff);

				if      (abs(flavId)==5) btsfutil.modifyBTagsWithSF(hasCSVtag_BTagDown, sfb - sfbunc,   beff);
				else if (abs(flavId)==4) btsfutil.modifyBTagsWithSF(hasCSVtag_BTagDown, sfb/5 - 2*sfbunc, beff);
				else                     btsfutil.modifyBTagsWithSF(hasCSVtag_BTagDown, sfl - sfbunc,   leff);
				}
			*/

			if(hasCSVtag || hasCSVtag_BTagUp || hasCSVtag_BTagDown)
				{
				selBJets.push_back(jet);
				}
			}

		if(debug){
			cout << "processed b-tagged jets" << endl;
			}


		/* Done before
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
		*/

		// -------------------------------------------------- all particles are selected

		if(debug){
			cout << "all particle-objects are processed, checking channel selection" << endl;
			}

		unsigned int n_leptons = selLeptons.size();
		// unsigned int n_taus = selTaus.size();
		unsigned int n_taus = selTausNoLep.size();
		//unsigned int n_jets = selJets.size();
		//unsigned int n_bjets = selBJets.size();
		// unsigned int n_jets = selSingleLepJets.size();
		unsigned int n_jets = selJetsNoLepNoTau.size();
		// unsigned int n_bjets = selSingleLepBJets.size();
		unsigned int n_bjets = selBJets.size();

		n_selected_leptons_weighted[n_leptons > 100 ? 99 : n_leptons] += weight;
		n_selected_taus_weighted [n_taus > 100 ? 99 : n_taus] += weight;
		n_selected_jets_weighted [n_jets > 100 ? 99 : n_jets] += weight;
		n_selected_bjets_weighted[n_bjets > 100 ? 99 : n_bjets] += weight;

		// and the sum of weight before splitting into channels:
		weight_before_channel_select += weight;
		increment(string("weight_before_channel_select"), weight);
		increment(string("weight_up_before_channel_select"), weight_up);
		increment(string("weight_down_before_channel_select"), weight_down);

		//
		// -------------------------------------------------- ASSIGN CHANNEL
		//

		// // TODO: did this classification above -- where the leptons are processed
		// // Event classification. Single lepton triggers are used for offline selection of dilepton events. The "else if"s guarantee orthogonality
		// if(selLeptons.size()>0)
		// 	slepId=selLeptons[0].pdgId();

		// bool 
		// 	isSingleMu(false),
		// 	isSingleE(false),
		// 	isDoubleMu(false),
		// 	isDoubleE(false),
		// 	isEMu(false);
		// // int multiChannel(0);


		// bool iso_lep = nVetoE==0 && nVetoMu==0 && selLeptons.size() == 1 && nGoodPV != 0; // 2^5
		// //if(selLeptons.size()!=1 || nGoodPV==0) continue; // Veto requirement alredy applied during the event categoriziation
		// isSingleMu = (abs(slepId)==13) && muTrigger && iso_lep;
		// isSingleE  = (abs(slepId)==11) && eTrigger  && iso_lep;
		// TODO: last discrepancy with multiselect!
		// TODO: and no double-lepton channel yet

		// --------------------------- store weights at different selections
		// Event selection booleans for el-tau and mu-tau channels


		// ------------------------------------------ Single lepton full analysis
		if (isSingleMu || isSingleE)
			{
			// in-channel selection/multiselect for leptons

			singlelep_ttbar_preselectedevents->Fill(1);


			// bool passJetRawSelection(selSingleLepJets.size()>1); // 2 jets
			//bool passJetSelection(selSingleLepJets.size()>1); // 2 jets // 2^4
			//bool passJetSelection(selJets.size()>1); // 2 jets // 2^4
			bool passJetSelection(n_jets>1); // 2 jets // 2^4
			bool passMetSelection(met.pt()>40.); // MET > 40 // 2^3
			// bool passMetSelection(n_met.pt()>40.); // MET > 40 // 2^3
			//bool passBtagsSelection(selSingleLepBJets.size()>0); // 1 b jet // 2^2
			//bool passBtagsSelection(selBJets.size()>0); // 1 b jet // 2^2
			bool passBtagsSelection(n_bjets>0); // 1 b jet // 2^2
			bool passTauSelection(n_taus==1); // only 1 tau // 2^1
			// bool passTauSelection(n_taus>0); // >= 1 tau in 13.2
			bool passOS( n_taus>0 && n_leptons>0 ? selLeptons[0].pdgId() * selTausNoLep[0].pdgId() < 0 : 0); // Oposite sign // 2^0


			// MULTISELECT
			// multiselection
			// TODO: multisel should be done per-channel, now it is one (el/mu) for all 
			unsigned int multisel = 0;
			// multisel += (isSingleMu ? 1 : 0); //! should be 1
			// multisel += (isSingleE ? 2 : 0);
			multisel += (passJetSelection ? 1 : 0);
			multisel += (passMetSelection ? 2 : 0);
			multisel += (passBtagsSelection ? 4 : 0);
			multisel += (passTauSelection ? 8 : 0);
			multisel += (passOS ? 16 : 0);

			fill_pt_e( string("singlelep_channel_met_pt"), met.pt(), weight);
			if (passJetSelection)
				{
				fill_pt_e( string("singlelep_jetsel_met_pt"), met.pt(), weight);
				}
			if (passJetSelection && passBtagsSelection && passTauSelection && passOS)
				{
				fill_pt_e( string("singlelep_allbutmetsel_met_pt"), met.pt(), weight);
				}


			if (passJetSelection && passMetSelection && passBtagsSelection && passTauSelection && passOS) {
				increment( string("weightflow_weight_passed_singlelep_selection"), weight );
				increment( string("weightflow_weight_up_passed_singlelep_selection"), weight_up );
				increment( string("weightflow_weight_down_passed_singlelep_selection"), weight_down );
				}

			if (isSingleMu)
				{
				fill_n( string("n_jets_singlemu"), n_jets, weight);
				fill_n( string("n_bjets_singlemu"), n_bjets, weight);
				fill_n( string("n_taus_singlemu"), n_taus, weight);

				fill_pt_e( string("singlemu_channel_met_pt"), met.pt(), weight);
				if (passJetSelection)
					{
					fill_pt_e( string("singlemu_jetsel_met_pt"), met.pt(), weight);
					}
				if (passJetSelection && passBtagsSelection && passTauSelection && passOS)
					{
					fill_pt_e( string("singlemu_allbutmetsel_met_pt"), met.pt(), weight);
					}

				weights_in_mu_channel[multisel] += weight;
				increment(string("weightflow_mu_") + to_string(multisel), weight);
				// increment(string("weightflow_up_mu_") + to_string(multisel), weight_up);
				// increment(string("weightflow_down_mu_") + to_string(multisel), weight_down);
				fill_pt_e( string("top1pt_muons_pt_individual"), selLeptons[0].pt(), weight);
				fill_pt_e( string("top1pt_muons_pt_individual_up"), selLeptons[0].pt(), weight_up);
				fill_pt_e( string("top1pt_muons_pt_individual_down"), selLeptons[0].pt(), weight_down);

				fill_particle_ids(string("nearest_particle_around_tau_singlemu"), closest_totaunolep_particle_id, weight);
				if (fabs(closest_totaunolep_particle_id) == 13)
					{
					increment(string("number_of_tausnolep_from_muon_found_in_singlemu"), 1);
					}

				if (passJetSelection)
					{
					fill_particle_ids(string("nearest_particle_around_tau_singlemu_jetselection"), closest_totaunolep_particle_id, weight);
					if (fabs(closest_totaunolep_particle_id) == 13)
						{
						increment(string("number_of_tausnolep_from_muon_found_in_singlemu_jetselection"), 1);
						}
					}

				if (passJetSelection && passMetSelection)
					{
					fill_n( string("singlemu_prebselpoint_n_jets"), n_jets, weight);
					fill_n( string("singlemu_prebselpoint_n_bjets"), n_bjets, weight);
					fill_n( string("singlemu_prebselpoint_n_taus"), n_taus, weight);
					}

				if (passJetSelection && passMetSelection && passBtagsSelection && passTauSelection && passOS)
					{
					fill_particle_ids(string("nearest_particle_around_tau_singlemu_fullselection"), closest_totaunolep_particle_id, weight);
					if (fabs(closest_totaunolep_particle_id) == 13)
						{
						increment(string("number_of_tausnolep_from_muon_found_in_singlemu_fullselection"), 1);
						}

					increment( string("weightflow_weight_passed_singlemu_selection"), weight );
					increment( string("weightflow_weight_up_passed_singlemu_selection"), weight_up );
					increment( string("weightflow_weight_down_passed_singlemu_selection"), weight_down );

					fill_n( string("singlemu_selection_n_jets"), n_jets, weight);
					fill_n( string("singlemu_selection_n_bjets"), n_bjets, weight);
					fill_n( string("singlemu_selection_n_taus"), n_taus, weight);

					// pts
					fill_pt_e( string("singlemu_selection_muon_pt"), selLeptons[0].pt(), weight);
					fill_pt_e( string("singlemu_selection_tau_pt"), selTausNoLep[0].pt(), weight);
					fill_pt_e( string("singlemu_selection_jet1_pt"), selJetsNoLepNoTau[0].pt(), weight);
					fill_pt_e( string("singlemu_selection_jet2_pt"), selJetsNoLepNoTau[1].pt(), weight);
					fill_pt_e( string("singlemu_selection_bjet_pt"), selBJets[0].pt(), weight);
					// fill_pt_e( string("singlemu_selection_met_pt"), n_met.pt(), weight);
					fill_pt_e( string("singlemu_selection_met_pt"), met.pt(), weight);

					// energies
					fill_pt_e( string("singlemu_selection_muon_energy"), selLeptons[0].energy(), weight);
					fill_pt_e( string("singlemu_selection_tau_energy"), selTausNoLep[0].energy(), weight);
					fill_pt_e( string("singlemu_selection_jet1_energy"), selJetsNoLepNoTau[0].energy(), weight);
					fill_pt_e( string("singlemu_selection_jet2_energy"), selJetsNoLepNoTau[1].energy(), weight);
					fill_pt_e( string("singlemu_selection_bjet_energy"), selBJets[0].energy(), weight);
					// fill_pt_e( string("singlemu_selection_met_energy"), n_met.energy(), weight);
					fill_pt_e( string("singlemu_selection_met_energy"), met.energy(), weight);

					// etas:
					fill_eta( string("singlemu_selection_muon_eta"), selLeptons[0].eta(), weight);
					fill_eta( string("singlemu_selection_tau_eta"), selTausNoLep[0].eta(), weight);
					fill_eta( string("singlemu_selection_jet1_eta"), selJetsNoLepNoTau[0].eta(), weight);
					fill_eta( string("singlemu_selection_jet2_eta"), selJetsNoLepNoTau[1].eta(), weight);
					fill_eta( string("singlemu_selection_bjet_eta"), selBJets[0].eta(), weight);
					// fill_eta( string("singlemu_selection_met_eta"), n_met.eta(), weight);
					fill_eta( string("singlemu_selection_met_eta"), met.eta(), weight);

					fill_pu( string("pileup_muselection_rawweight_pergoodpv"), nGoodPV, rawWeight);
					fill_pu( string("pileup_muselection_weight_pergoodpv"), nGoodPV, weight);
					fill_pu( string("pileup_muselection_weight_up_pergoodpv"), nGoodPV, weight_up);
					fill_pu( string("pileup_muselection_weight_down_pergoodpv"), nGoodPV, weight_down);

					fill_pu( string("pileup_muselection_rawweight_pernuminters"), num_inters, rawWeight);
					fill_pu( string("pileup_muselection_weight_pernuminters"), num_inters, weight);
					fill_pu( string("pileup_muselection_weight_up_pernuminters"), num_inters, weight_up);
					fill_pu( string("pileup_muselection_weight_down_pernuminters"), num_inters, weight_down);
					}
				}

			if (isSingleE)
				{
				fill_n( string("n_jets_singleel"), n_jets, weight);
				fill_n( string("n_bjets_singleel"), n_bjets, weight);
				fill_n( string("n_taus_singleel"), n_taus, weight);

				fill_pt_e( string("singleel_channel_met_pt"), met.pt(), weight);
				if (passJetSelection)
					{
					fill_pt_e( string("singleel_jetsel_met_pt"), met.pt(), weight);
					}
				if (passJetSelection && passBtagsSelection && passTauSelection && passOS)
					{
					fill_pt_e( string("singleel_allbutmetsel_met_pt"), met.pt(), weight);
					}

				weights_in_el_channel[multisel] += weight;
				increment(string("weightflow_e_") + to_string(multisel), weight);
				// increment(string("weightflow_up_e_") + to_string(multisel), weight_up);
				// increment(string("weightflow_down_e_") + to_string(multisel), weight_down);
				fill_pt_e( string("top1pt_electrons_pt_individual"), selLeptons[0].pt(), weight);
				fill_pt_e( string("top1pt_electrons_pt_individual_up"), selLeptons[0].pt(), weight_up);
				fill_pt_e( string("top1pt_electrons_pt_individual_down"), selLeptons[0].pt(), weight_down);

				fill_particle_ids(string("nearest_particle_around_tau_singleel"), closest_totaunolep_particle_id, weight);
				if (fabs(closest_totaunolep_particle_id) == 11)
					{
					increment(string("number_of_tausnolep_from_electron_found_in_singleel"), 1);
					}

				if (passJetSelection)
					{
					fill_particle_ids(string("nearest_particle_around_tau_singleel_jetselection"), closest_totaunolep_particle_id, weight);
					if (fabs(closest_totaunolep_particle_id) == 11)
						{
						increment(string("number_of_tausnolep_from_electron_found_in_singleel_jetselection"), 1);
						}
					}

				if (passJetSelection && passMetSelection)
					{
					fill_n( string("singleel_prebselpoint_n_jets"), n_jets, weight);
					fill_n( string("singleel_prebselpoint_n_bjets"), n_bjets, weight);
					fill_n( string("singleel_prebselpoint_n_taus"), n_taus, weight);
					}

				if (passJetSelection && passMetSelection && passBtagsSelection && passTauSelection && passOS)
					{
					fill_particle_ids(string("nearest_particle_around_tau_singleel_fullselection"), closest_totaunolep_particle_id, weight);
					if (fabs(closest_totaunolep_particle_id) == 11)
						{
						increment(string("number_of_tausnolep_from_electron_found_in_singleel_fullselection"), 1);
						}

					increment( string("weightflow_weight_passed_singleel_selection"), weight );
					increment( string("weightflow_weight_up_passed_singleel_selection"), weight_up );
					increment( string("weightflow_weight_down_passed_singleel_selection"), weight_down );

					fill_n( string("singleel_selection_n_jets"), n_jets, weight);
					fill_n( string("singleel_selection_n_bjets"), n_bjets, weight);
					fill_n( string("singleel_selection_n_taus"), n_taus, weight);

					fill_pt_e( string("singleel_selection_electron_pt"), selLeptons[0].pt(), weight);
					fill_pt_e( string("singleel_selection_tau_pt"), selTausNoLep[0].pt(), weight);
					fill_pt_e( string("singleel_selection_jet1_pt"), selJetsNoLepNoTau[0].pt(), weight);
					fill_pt_e( string("singleel_selection_jet2_pt"), selJetsNoLepNoTau[1].pt(), weight);
					fill_pt_e( string("singleel_selection_bjet_pt"), selBJets[0].pt(), weight);
					// fill_pt_e( string("singleel_selection_met_pt"), n_met.pt(), weight);
					fill_pt_e( string("singlemu_selection_met_pt"), met.pt(), weight);

					// energies
					fill_pt_e( string("singleel_selection_electron_energy"), selLeptons[0].energy(), weight);
					fill_pt_e( string("singleel_selection_tau_energy"), selTausNoLep[0].energy(), weight);
					fill_pt_e( string("singleel_selection_jet1_energy"), selJetsNoLepNoTau[0].energy(), weight);
					fill_pt_e( string("singleel_selection_jet2_energy"), selJetsNoLepNoTau[1].energy(), weight);
					fill_pt_e( string("singleel_selection_bjet_energy"), selBJets[0].energy(), weight);
					// fill_pt_e( string("singleel_selection_met_energy"), n_met.energy(), weight);
					fill_pt_e( string("singlemu_selection_met_energy"), met.energy(), weight);

					//etas:
					fill_eta( string("singleel_selection_electron_eta"), selLeptons[0].eta(), weight);
					fill_eta( string("singleel_selection_tau_eta"), selTausNoLep[0].eta(), weight);
					fill_eta( string("singleel_selection_jet1_eta"), selJetsNoLepNoTau[0].eta(), weight);
					fill_eta( string("singleel_selection_jet2_eta"), selJetsNoLepNoTau[1].eta(), weight);
					fill_eta( string("singleel_selection_bjet_eta"), selBJets[0].eta(), weight);
					// fill_eta( string("singleel_selection_met_eta"), n_met.eta(), weight);
					fill_eta( string("singlemu_selection_met_eta"), met.eta(), weight);

					fill_pu( string("pileup_elselection_rawweight_pergoodpv"), nGoodPV, rawWeight);
					fill_pu( string("pileup_elselection_weight_pergoodpv"), nGoodPV, weight);
					fill_pu( string("pileup_elselection_weight_up_pergoodpv"), nGoodPV, weight_up);
					fill_pu( string("pileup_elselection_weight_down_pergoodpv"), nGoodPV, weight_down);

					fill_pu( string("pileup_elselection_rawweight_pernuminters"), num_inters, rawWeight);
					fill_pu( string("pileup_elselection_weight_pernuminters"), num_inters, weight);
					fill_pu( string("pileup_elselection_weight_up_pernuminters"), num_inters, weight_up);
					fill_pu( string("pileup_elselection_weight_down_pernuminters"), num_inters, weight_down);
					}
				}

			// if(passJetSelection && passBtagsSelection) // 2 jets, 1 b jet, 1 isolated lepton
			// 	{
			// 	 now these histograms are dissabled -- use counters to substitute them
			// 	if(isSingleMu) singlelep_ttbar_selected_mu_events->Fill(1);
			// 	else if (isSingleE) singlelep_ttbar_selected_el_events->Fill(1);
			// 	crossel_sum_weights_raw += rawWeight;
			// 	crossel_sum_weights += weight;
			// 	}

			/*
			if(passJetSelection && passMetSelection && passBtagsSelection && passTauSelection && passOS )
				{
				// TODO: try saving the whole event
				if(isSingleMu)
					{
					//singlelep_ttbar_selected_mu_events->Fill(1);
					oursel_sum_weights_mu += weight;
					}
				else if (isSingleE)
					{
					//singlelep_ttbar_selected_el_events->Fill(1);
					oursel_sum_weights_el += weight;
					}
				fprintf(csv_out, "oursel:%d,%d,%d,%g,%g,%g,%g,%d,", iev, num_inters, nGoodPV, rawWeight, weight, weight_up, weight_down, isSingleE);
				oursel_sum_weights_raw += rawWeight;
				oursel_sum_weights += weight;

				// METs with corrections
				// LorentzVector met_values[7];
				fprintf(csv_out, "%g,%g,%g,%g,%g,%g,%g,", met_pt_values[0], met_pt_values[1], met_pt_values[2], met_pt_values[3], met_pt_values[4], met_pt_values[5], met_pt_values[6]);
				fprintf(csv_out, "%g,", selLeptons[0].vz());
				fprintf(csv_out, "%g,%g,%g,%g,",  selLeptons[0].px(), selLeptons[0].py(), selLeptons[0].pz(), selLeptons[0].pt());
				//fprintf(csv_out, "%g,", selTausNoLep[0].vz());
				// selTausNoLep
				fprintf(csv_out, "%g,", selTausNoLep[0].vz());
				fprintf(csv_out, "%g,%g,%g,%g,", selTausNoLep[0].px(), selTausNoLep[0].py(), selTausNoLep[0].pz(), selTausNoLep[0].pt() );

				// selBJets
				fprintf(csv_out, "%g,", selBJets[0].vz());
				fprintf(csv_out, "%g,%g,%g,%g,",  selBJets[0].px(), selBJets[0].py(), selBJets[0].pz(), selBJets[0].pt());

				// selJetsNoLepNoTau
				//fprintf(csv_out, "%g,", selJetsNoLepNoTau[0].vz());
				fprintf(csv_out, "%g,", selJetsNoLepNoTau[0].vz());
				fprintf(csv_out, "%g,%g,%g,%g,", selJetsNoLepNoTau[0].px(), selJetsNoLepNoTau[0].py(), selJetsNoLepNoTau[0].pz(), selJetsNoLepNoTau[0].pt() );
				//fprintf(csv_out, "%g,", selJetsNoLepNoTau[1].vz());
				fprintf(csv_out, "%g,", selJetsNoLepNoTau[1].vz());
				fprintf(csv_out, "%g,%g,%g,%g\n", selJetsNoLepNoTau[1].px(), selJetsNoLepNoTau[1].py(), selJetsNoLepNoTau[1].pz(), selJetsNoLepNoTau[1].pt() );

				}
				*/
			}

		// inline control functions usage:
		//   fill_pt_e( "control_point_name", value, weight)
		// FIXME: do NEWMULTISELECT somehow well
		//   fill_eta( "control_point_name", value, weight )   <-- different TH1F range and binning
		//   increment( "control_point_name", weight )


		// --------------------------------------------- Dilepton full analysis
		if (isDoubleE || isDoubleMu || isEMu)
			{
			int dilId (1);
			// slepId(0);
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

			// Event selection booleans
			bool passMllVeto(isEMu ? dileptonSystem.mass()>12. : (fabs(dileptonSystem.mass()-91.)>15 && dileptonSystem.mass()>12. ) );
			// bool passJetSelection(selJets.size()>1);
			bool passJetSelection(n_jets>1);
			bool passMetSelection(met.pt()>40.);
			bool passOS(selLeptons[0].pdgId() * selLeptons[1].pdgId() < 0 );
			// bool passBtagsSelection(selBJets.size()>1); // FIXME: differentiate chhiggs selection from cross-section selection
			bool passBtagsSelection(n_bjets>0);

			// MULTISELECT
			unsigned int multisel = 0;
			multisel += (passMllVeto ? 1 : 0);
			multisel += (passJetSelection ? 2 : 0);
			multisel += (passMetSelection ? 4 : 0);
			multisel += (passOS ? 8 : 0);
			multisel += (passBtagsSelection ? 16 : 0);

			fill_pt_e( string("doublelep_channel_met_pt"), met.pt(), weight);
			if (passMllVeto && passJetSelection)
				{
				fill_pt_e( string("doublelep_jetsel_met_pt"), met.pt(), weight);
				}
			if (passMllVeto && passJetSelection && passOS && passBtagsSelection)
				{
				fill_pt_e( string("doublelep_allbutmetsel_met_pt"), met.pt(), weight);
				}

			if (passMllVeto && passJetSelection && passMetSelection && passOS && passBtagsSelection)
				{
				increment( string("weightflow_weight_passed_doublelep_selection"), weight );
				increment( string("weightflow_weight_up_passed_doublelep_selection"), weight_up );
				increment( string("weightflow_weight_down_passed_doublelep_selection"), weight_down );
				}

			if (isDoubleE)
				{
				increment(string("weightflow_ee_") + to_string(multisel), weight);
				// increment(string("weightflow_up_ee_") + to_string(multisel), weight_up);
				// increment(string("weightflow_down_ee_") + to_string(multisel), weight_down);

				fill_pt_e( string("doubleel_channel_met_pt"), met.pt(), weight);
				if (passMllVeto && passJetSelection)
					{
					fill_pt_e( string("doubleel_jetsel_met_pt"), met.pt(), weight);
					}
				if (passMllVeto && passJetSelection && passOS && passBtagsSelection)
					{
					fill_pt_e( string("doubleel_allbutmetsel_met_pt"), met.pt(), weight);
					}

				fill_n( string("elel_channel_n_jets"), n_jets, weight);
				fill_n( string("elel_channel_n_bjets"), n_bjets, weight);
				fill_n( string("elel_channel_n_taus"), n_taus, weight);

				if( passMllVeto && passJetSelection && passMetSelection && passOS)
					{
					fill_n( string("elel_prebselpoint_n_jets"), n_jets, weight);
					fill_n( string("elel_prebselpoint_n_bjets"), n_bjets, weight);
					fill_n( string("elel_prebselpoint_n_taus"), n_taus, weight);
					}

				if( passMllVeto && passJetSelection && passMetSelection && passOS && passBtagsSelection)
					{
					increment( string("weightflow_weight_passed_doubleel_selection"), weight );

					fill_n( string("elel_selection_n_jets"), n_jets, weight);
					fill_n( string("elel_selection_n_bjets"), n_bjets, weight);
					fill_n( string("elel_selection_n_taus"), n_taus, weight);

					fill_pt_e( string("elel_selection_el1_pt"), selElectrons[0].pt(), weight);
					fill_pt_e( string("elel_selection_el2_pt"), selElectrons[1].pt(), weight);
					// fill_pt_e( string("elel_selection_tau_pt"), selTausNoLep[0].pt(), weight);
					fill_pt_e( string("elel_selection_jet1_pt"), selJetsNoLepNoTau[0].pt(), weight);
					fill_pt_e( string("elel_selection_jet2_pt"), selJetsNoLepNoTau[1].pt(), weight);
					fill_pt_e( string("elel_selection_bjet_pt"), selBJets[0].pt(), weight);
					// fill_pt_e( string("elel_selection_met_pt"), n_met.pt(), weight);
					fill_pt_e( string("elel_selection_met_pt"), met.pt(), weight);

					// energies
					fill_pt_e( string("elel_selection_el1_energy"), selElectrons[0].energy(), weight);
					fill_pt_e( string("elel_selection_el2_energy"), selElectrons[1].energy(), weight);
					// fill_pt_e( string("elel_selection_tau_energy"), selTausNoLep[0].energy(), weight);
					fill_pt_e( string("elel_selection_jet1_energy"), selJetsNoLepNoTau[0].energy(), weight);
					fill_pt_e( string("elel_selection_jet2_energy"), selJetsNoLepNoTau[1].energy(), weight);
					fill_pt_e( string("elel_selection_bjet_energy"), selBJets[0].energy(), weight);
					// fill_pt_e( string("elel_selection_met_energy"), n_met.energy(), weight);
					fill_pt_e( string("elel_selection_met_energy"), met.energy(), weight);

					//etas:
					fill_eta( string("elel_selection_el1_eta"), selElectrons[0].eta(), weight);
					fill_eta( string("elel_selection_el2_eta"), selElectrons[1].eta(), weight);
					// fill_eta( string("elel_selection_tau_eta"), selTausNoLep[0].eta(), weight);
					fill_eta( string("elel_selection_jet1_eta"), selJetsNoLepNoTau[0].eta(), weight);
					fill_eta( string("elel_selection_jet2_eta"), selJetsNoLepNoTau[1].eta(), weight);
					fill_eta( string("elel_selection_bjet_eta"), selBJets[0].eta(), weight);
					// fill_eta( string("elel_selection_met_eta"), n_met.eta(), weight);
					fill_eta( string("elel_selection_met_eta"), met.eta(), weight);

					fill_pu( string("pileup_elelselection_rawweight_pergoodpv"), nGoodPV, rawWeight);
					fill_pu( string("pileup_elelselection_weight_pergoodpv"), nGoodPV, weight);
					fill_pu( string("pileup_elelselection_weight_up_pergoodpv"), nGoodPV, weight_up);
					fill_pu( string("pileup_elelselection_weight_down_pergoodpv"), nGoodPV, weight_down);

					fill_pu( string("pileup_elelselection_rawweight_pernuminters"), num_inters, rawWeight);
					fill_pu( string("pileup_elelselection_weight_pernuminters"), num_inters, weight);
					fill_pu( string("pileup_elelselection_weight_up_pernuminters"), num_inters, weight_up);
					fill_pu( string("pileup_elelselection_weight_down_pernuminters"), num_inters, weight_down);
					}
				}

			if (isDoubleMu)
				{
				increment(string("weightflow_mumu_") + to_string(multisel), weight);
				// increment(string("weightflow_up_mumu_") + to_string(multisel), weight_up);
				// increment(string("weightflow_down_mumu_") + to_string(multisel), weight_down);

				fill_pt_e( string("doublemu_channel_met_pt"), met.pt(), weight);
				if (passMllVeto && passJetSelection)
					{
					fill_pt_e( string("doublemu_jetsel_met_pt"), met.pt(), weight);
					}
				if (passMllVeto && passJetSelection && passOS && passBtagsSelection)
					{
					fill_pt_e( string("doublemu_allbutmetsel_met_pt"), met.pt(), weight);
					}

				fill_n( string("mumu_channel_n_jets"), n_jets, weight);
				fill_n( string("mumu_channel_n_bjets"), n_bjets, weight);
				fill_n( string("mumu_channel_n_taus"), n_taus, weight);

				if(passMllVeto && passJetSelection && passMetSelection && passOS)
					{
					fill_n( string("mumu_prebselpoint_n_jets"), n_jets, weight);
					fill_n( string("mumu_prebselpoint_n_bjets"), n_bjets, weight);
					fill_n( string("mumu_prebselpoint_n_taus"), n_taus, weight);
					}

				if(passMllVeto && passJetSelection && passMetSelection && passOS && passBtagsSelection)
					{
					increment( string("weightflow_weight_passed_doublemu_selection"), weight );

					fill_n( string("mumu_selection_n_jets"), n_jets, weight);
					fill_n( string("mumu_selection_n_bjets"), n_bjets, weight);
					fill_n( string("mumu_selection_n_taus"), n_taus, weight);

					fill_pt_e( string("mumu_selection_mu1_pt"), selMuons[0].pt(), weight);
					fill_pt_e( string("mumu_selection_mu2_pt"), selMuons[1].pt(), weight);
					// fill_pt_e( string("mumu_selection_tau_pt"), selTausNoLep[0].pt(), weight);
					fill_pt_e( string("mumu_selection_jet1_pt"), selJetsNoLepNoTau[0].pt(), weight);
					fill_pt_e( string("mumu_selection_jet2_pt"), selJetsNoLepNoTau[1].pt(), weight);
					fill_pt_e( string("mumu_selection_bjet_pt"), selBJets[0].pt(), weight);
					// fill_pt_e( string("mumu_selection_met_pt"), n_met.pt(), weight);
					fill_pt_e( string("mumu_selection_met_pt"), met.pt(), weight);

					// energies
					fill_pt_e( string("mumu_selection_mu1_energy"), selMuons[0].energy(), weight);
					fill_pt_e( string("mumu_selection_mu2_energy"), selMuons[1].energy(), weight);
					// fill_pt_e( string("mumu_selection_tau_energy"), selTausNoLep[0].energy(), weight);
					fill_pt_e( string("mumu_selection_jet1_energy"), selJetsNoLepNoTau[0].energy(), weight);
					fill_pt_e( string("mumu_selection_jet2_energy"), selJetsNoLepNoTau[1].energy(), weight);
					fill_pt_e( string("mumu_selection_bjet_energy"), selBJets[0].energy(), weight);
					// fill_pt_e( string("mumu_selection_met_energy"), n_met.energy(), weight);
					fill_pt_e( string("mumu_selection_met_energy"), met.energy(), weight);

					// etas
					fill_eta( string("mumu_selection_mu1_eta"), selMuons[0].eta(), weight);
					fill_eta( string("mumu_selection_mu2_eta"), selMuons[1].eta(), weight);
					// fill_eta( string("mumu_selection_tau_eta"), selTausNoLep[0].eta(), weight);
					fill_eta( string("mumu_selection_jet1_eta"), selJetsNoLepNoTau[0].eta(), weight);
					fill_eta( string("mumu_selection_jet2_eta"), selJetsNoLepNoTau[1].eta(), weight);
					fill_eta( string("mumu_selection_bjet_eta"), selBJets[0].eta(), weight);
					// fill_eta( string("mumu_selection_met_eta"), n_met.eta(), weight);
					fill_eta( string("mumu_selection_met_eta"), met.eta(), weight);

					fill_pu( string("pileup_mumuselection_rawweight_pergoodpv"), nGoodPV, rawWeight);
					fill_pu( string("pileup_mumuselection_weight_pergoodpv"), nGoodPV, weight);
					fill_pu( string("pileup_mumuselection_weight_up_pergoodpv"), nGoodPV, weight_up);
					fill_pu( string("pileup_mumuselection_weight_down_pergoodpv"), nGoodPV, weight_down);

					fill_pu( string("pileup_mumuselection_rawweight_pernuminters"), num_inters, rawWeight);
					fill_pu( string("pileup_mumuselection_weight_pernuminters"), num_inters, weight);
					fill_pu( string("pileup_mumuselection_weight_up_pernuminters"), num_inters, weight_up);
					fill_pu( string("pileup_mumuselection_weight_down_pernuminters"), num_inters, weight_down);
					}
				}

			if (isEMu)
				{
				increment(string("weightflow_emu_") + to_string(multisel), weight);
				// increment(string("weightflow_up_emu_") + to_string(multisel), weight_up);
				// increment(string("weightflow_down_emu_") + to_string(multisel), weight_down);

				fill_pt_e( string("doubleelmu_channel_met_pt"), met.pt(), weight);
				if (passMllVeto && passJetSelection)
					{
					fill_pt_e( string("doubleelmu_jetsel_met_pt"), met.pt(), weight);
					}
				if (passMllVeto && passJetSelection && passOS && passBtagsSelection)
					{
					fill_pt_e( string("doubleelmu_allbutmetsel_met_pt"), met.pt(), weight);
					}

				fill_n( string("elmu_channel_n_jets"), n_jets, weight);
				fill_n( string("elmu_channel_n_bjets"), n_bjets, weight);
				fill_n( string("elmu_channel_n_taus"), n_taus, weight);

				if (passMllVeto && passJetSelection && passMetSelection && passOS)
					{
					fill_n( string("elmu_prebselpoint_n_jets"), n_jets, weight);
					fill_n( string("elmu_prebselpoint_n_bjets"), n_bjets, weight);
					fill_n( string("elmu_prebselpoint_n_taus"), n_taus, weight);
					}

				if (passMllVeto && passJetSelection && passMetSelection && passOS && passBtagsSelection)
					{
					increment( string("weightflow_weight_passed_doubleemu_selection"), weight );

					fill_n( string("elmu_selection_n_jets"), n_jets, weight);
					fill_n( string("elmu_selection_n_bjets"), n_bjets, weight);
					fill_n( string("elmu_selection_n_taus"), n_taus, weight);

					fill_pt_e( string("elmu_selection_el_pt"), selElectrons[0].pt(), weight);
					fill_pt_e( string("elmu_selection_mu_pt"), selMuons[0].pt(), weight);
					// fill_pt_e( string("elmu_selection_tau_pt"), selTausNoLep[0].pt(), weight);
					fill_pt_e( string("elmu_selection_jet1_pt"), selJetsNoLepNoTau[0].pt(), weight);
					fill_pt_e( string("elmu_selection_jet2_pt"), selJetsNoLepNoTau[1].pt(), weight);
					fill_pt_e( string("elmu_selection_bjet_pt"), selBJets[0].pt(), weight);
					// fill_pt_e( string("elmu_selection_met_pt"), n_met.pt(), weight);
					fill_pt_e( string("elmu_selection_met_pt"), met.pt(), weight);

					// energies
					fill_pt_e( string("elmu_selection_el_energy"), selElectrons[0].energy(), weight);
					fill_pt_e( string("elmu_selection_mu_energy"), selMuons[0].energy(), weight);
					// fill_pt_e( string("elmu_selection_tau_energy"), selTausNoLep[0].energy(), weight);
					fill_pt_e( string("elmu_selection_jet1_energy"), selJetsNoLepNoTau[0].energy(), weight);
					fill_pt_e( string("elmu_selection_jet2_energy"), selJetsNoLepNoTau[1].energy(), weight);
					fill_pt_e( string("elmu_selection_bjet_energy"), selBJets[0].energy(), weight);
					// fill_pt_e( string("elmu_selection_met_energy"), n_met.energy(), weight);
					fill_pt_e( string("elmu_selection_met_energy"), met.energy(), weight);

					// etas
					fill_eta( string("elmu_selection_el_eta"), selElectrons[0].eta(), weight);
					fill_eta( string("elmu_selection_mu_eta"), selMuons[0].eta(), weight);
					// fill_eta( string("elmu_selection_tau_eta"), selTausNoLep[0].eta(), weight);
					fill_eta( string("elmu_selection_jet1_eta"), selJetsNoLepNoTau[0].eta(), weight);
					fill_eta( string("elmu_selection_jet2_eta"), selJetsNoLepNoTau[1].eta(), weight);
					fill_eta( string("elmu_selection_bjet_eta"), selBJets[0].eta(), weight);
					// fill_eta( string("elmu_selection_met_eta"), n_met.eta(), weight);
					fill_eta( string("elmu_selection_met_eta"), met.eta(), weight);

					fill_pu( string("pileup_elmuselection_rawweight_pergoodpv"), nGoodPV, rawWeight);
					fill_pu( string("pileup_elmuselection_weight_pergoodpv"), nGoodPV, weight);
					fill_pu( string("pileup_elmuselection_weight_up_pergoodpv"), nGoodPV, weight_up);
					fill_pu( string("pileup_elmuselection_weight_down_pergoodpv"), nGoodPV, weight_down);

					fill_pu( string("pileup_elmuselection_rawweight_pernuminters"), num_inters, rawWeight);
					fill_pu( string("pileup_elmuselection_weight_pernuminters"), num_inters, weight);
					fill_pu( string("pileup_elmuselection_weight_up_pernuminters"), num_inters, weight_up);
					fill_pu( string("pileup_elmuselection_weight_down_pernuminters"), num_inters, weight_down);
					}
				}
			}

		/* debugging
		if (multisel > MULTISEL_SIZE)
			{
			printf("in event %d too large multisel: %d", iev, multisel);
			break;
			}
		*/

		// TODO: properly count multichannel?
		if (isSingleE && isSingleMu) nMultiChannel++;

		if(debug){
			cout << "channel is defined, running the event selection" << endl;
			}


		} // End single file event loop

	fprintf(csv_out, "acceptances:");
	fprintf(csv_out, "%s,", urls[f].c_str());
	fprintf(csv_out, "%d,%d,%g,%g,", iev, n_events_pass_lumi, sum_weights_raw, sum_weights);
	fprintf(csv_out, "%g,%g,", sum_weights_passtrig_raw, sum_weights_passtrig);
	fprintf(csv_out, "%g,%g,",  crossel_sum_weights_raw, crossel_sum_weights);
	fprintf(csv_out, "%g,%g,",  oursel_sum_weights_raw, oursel_sum_weights);
	fprintf(csv_out, "%g,%g,",  oursel_sum_weights_el, oursel_sum_weights_mu);
	fprintf(csv_out, "%g,%g\n", marasel_sum_weights_raw, marasel_sum_weights);

	// MULTISELECT

	fprintf(csv_out, "weights_in_no_channel, %d, %d, %g, %g, %g, %g", iev, n_events_pass_lumi, sum_weights, sum_weights_passtrig_raw, sum_weights_passtrig, weight_before_channel_select);
	for (int i=0; i<MULTISEL_SIZE; i++)
		{
		fprintf(csv_out, ",%g", weights_in_no_channel[i]);
		//fprintf(csv_out, "%d,", weights_in_selections_int[i]);
		}
	fprintf(csv_out, "\n");

	fprintf(csv_out, "weights_in_el_channel, %d, %d, %g, %g, %g, %g", iev, n_events_pass_lumi, sum_weights, sum_weights_passtrig_raw, sum_weights_passtrig, weight_before_channel_select);
	for (int i=0; i<MULTISEL_SIZE; i++)
		{
		fprintf(csv_out, ",%g", weights_in_el_channel[i]);
		//fprintf(csv_out, "%d,", weights_in_selections_int[i]);
		}
	fprintf(csv_out, "\n");

	fprintf(csv_out, "weights_in_mu_channel, %d, %d, %g, %g, %g, %g", iev, n_events_pass_lumi, sum_weights, sum_weights_passtrig_raw, sum_weights_passtrig, weight_before_channel_select);
	for (int i=0; i<MULTISEL_SIZE; i++)
		{
		fprintf(csv_out, ",%g", weights_in_mu_channel[i]);
		//fprintf(csv_out, "%d,", weights_in_selections_int[i]);
		}
	fprintf(csv_out, "\n");

	fprintf(csv_out, "weights_in_elmu_channel, %d, %d, %g, %g, %g, %g", iev, n_events_pass_lumi, sum_weights, sum_weights_passtrig_raw, sum_weights_passtrig, weight_before_channel_select);
	for (int i=0; i<MULTISEL_SIZE; i++)
		{
		fprintf(csv_out, ",%g", weights_in_elmu_channel[i]);
		//fprintf(csv_out, "%d,", weights_in_selections_int[i]);
		}
	fprintf(csv_out, "\n");

	fprintf(csv_out, "weights_in_elel_channel, %d, %d, %g, %g, %g, %g", iev, n_events_pass_lumi, sum_weights, sum_weights_passtrig_raw, sum_weights_passtrig, weight_before_channel_select);
	for (int i=0; i<MULTISEL_SIZE; i++)
		{
		fprintf(csv_out, ",%g", weights_in_elel_channel[i]);
		//fprintf(csv_out, "%d,", weights_in_selections_int[i]);
		}
	fprintf(csv_out, "\n");

	fprintf(csv_out, "weights_in_mumu_channel, %d, %d, %g, %g, %g, %g", iev, n_events_pass_lumi, sum_weights, sum_weights_passtrig_raw, sum_weights_passtrig, weight_before_channel_select);
	for (int i=0; i<MULTISEL_SIZE; i++)
		{
		fprintf(csv_out, ",%g", weights_in_mumu_channel[i]);
		//fprintf(csv_out, "%d,", weights_in_selections_int[i]);
		}
	fprintf(csv_out, "\n");




	fprintf(csv_out, "negative_events_nvtx:");
	for (int i=0; i<100; i++)
		{ fprintf(csv_out, "%d,", negative_event_nvtx[i]); }
	fprintf(csv_out, "\n");
	fprintf(csv_out, "positive_events_nvtx:");
	for (int i=0; i<100; i++)
		{ fprintf(csv_out, "%d,", positive_event_nvtx[i]); }
	fprintf(csv_out, "\n");

	fprintf(csv_out, "negative_event_pernvtx_weight:");
	for (int i=0; i<100; i++)
		{ fprintf(csv_out, "%g,", negative_event_pernvtx_weight[i]); }
	fprintf(csv_out, "\n");

	fprintf(csv_out, "positive_event_pernvtx_weight:");
	for (int i=0; i<100; i++)
		{ fprintf(csv_out, "%g,", positive_event_pernvtx_weight[i]); }
	fprintf(csv_out, "\n");


	// double event_pergoodpv_weight[100];
	fprintf(csv_out, "event_pergoodpv_weight:");
	for (int i=0; i<100; i++)
		{ fprintf(csv_out, "%g,", event_pergoodpv_weight[i]); }
	fprintf(csv_out, "\n");

	// double negative_event_pergoodpv_weight[100];
	fprintf(csv_out, "negative_event_pergoodpv_weight:");
	for (int i=0; i<100; i++)
		{ fprintf(csv_out, "%g,", negative_event_pergoodpv_weight[i]); }
	fprintf(csv_out, "\n");

	// double positive_event_pergoodpv_weight[100];
	fprintf(csv_out, "positive_event_pergoodpv_weight:");
	for (int i=0; i<100; i++)
		{ fprintf(csv_out, "%g,", positive_event_pergoodpv_weight[i]); }
	fprintf(csv_out, "\n");

	// double n_selected_leptons_weighted[100];
	fprintf(csv_out, "n_selected_leptons_weighted:");
	for (int i=0; i<100; i++)
		{ fprintf(csv_out, "%g,", n_selected_leptons_weighted[i]); }
	fprintf(csv_out, "\n");

	// double n_selected_taus_weighted[100];
	fprintf(csv_out, "n_selected_taus_weighted:");
	for (int i=0; i<100; i++)
		{ fprintf(csv_out, "%g,", n_selected_taus_weighted[i]); }
	fprintf(csv_out, "\n");

	// double n_selected_jets_weighted[100];
	fprintf(csv_out, "n_selected_jets_weighted:");
	for (int i=0; i<100; i++)
		{ fprintf(csv_out, "%g,", n_selected_jets_weighted[i]); }
	fprintf(csv_out, "\n");

	// double n_selected_bjets_weighted[100];
	fprintf(csv_out, "n_selected_bjets_weighted:");
	for (int i=0; i<100; i++)
		{ fprintf(csv_out, "%g,", n_selected_bjets_weighted[i]); }
	fprintf(csv_out, "\n");

	printf("\n");
	printf("Done processing the file\n");
	printf("\n");

	delete file;

	fprintf(csv_out, "End of event loop in the file.\n\n");
	} // End loop on files

printf("Done processing the job of files\n");

fprintf(csv_out, "End of (file loop) the job.\n");

// Controls distributions of processed particles


// CONTROLINFO

fprintf(csv_out, "New output (sums per whole job!):\n");

// inline control functions usage:
//   fill_pt_e( "control_point_name", value, weight)
//   fill_eta( "control_point_name", value, weight )   <-- different TH1F range and binning
//   increment( "control_point_name", weight )
//   printout_distrs(FILE * out)
//   printout_counters(FILE * out)

//std::map<string, TH1D*> th1d_distr_control;
//top2pt_jets_pt_taucleaned

// cout << "some control distrs:" << endl;
// cout << th1d_distr_control["top2pt_jets_pt_taucleaned"].Integral() << endl;
// cout << th1d_distr_control["top2pt_jets_pt_taucleaned"].GetSize() << endl;

// cout << th1d_distr_control["all_jets_pt_taucleaned"].Integral() << endl;
// cout << th1d_distr_control["all_jets_pt_taucleaned"].GetSize() << endl;

// cout << th1d_distr_control["all_jets_pt_correctedF"].Integral() << endl;
// cout << th1d_distr_control["all_jets_pt_correctedF"].GetSize() << endl;

// cout << th1d_distr_control["all_jets_pt_corrected2"].Integral() << endl;
// cout << th1d_distr_control["all_jets_pt_corrected2"].GetSize() << endl;

// cout << th1d_distr_control["all_jets_pt_slimmed"].Integral() << endl;
// cout << th1d_distr_control["all_jets_pt_slimmed"].GetSize() << endl;

// cout << th1d_distr_control["all_jets_pt_corrected1"].Integral() << endl;
// cout << th1d_distr_control["all_jets_pt_corrected1"].GetSize() << endl;

// hopefully TString will get converted to string...

//printout_counters(csv_out, string(isMC ? "MC,": "Data,") + dtag_s + string(",") + job_num);
//printout_distrs(csv_out, string(isMC ? "MC,": "Data,") + dtag_s + string(",") + job_num);
printout_counters(csv_out, job_def);
printout_distrs(csv_out, job_def);

// So, each job output contains for each value:
// value_name,MC/Data,dtag,job_num,value
// and the table for distr:
// distr_name,MC/Data,dtag,job_num,bin,contents
// simple grep/ack by first coloumn produces a ready data.frame for R
// which should be summarised for each dataset -- easy with ddply
// the only problem is:
// FIXME: how to do MULTISELECT nicely here?


fprintf(csv_out, "End of job output.\n\n");


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

// re-enabled the ROOT output
// for resubmit option of the job script
TFile *ofile = TFile::Open (outUrl + ".root", "recreate");

singlelep_ttbar_initialevents->Write();
singlelep_ttbar_preselectedevents->Write();
// singlelep_ttbar_selected_mu_events->Write();
// singlelep_ttbar_selected_el_events->Write();
// singlelep_ttbar_selected2_mu_events->Write();
// singlelep_ttbar_selected2_el_events->Write();

// singlelep_ttbar_maraselected_mu_events->Write();
// singlelep_ttbar_maraselected_el_events->Write();

ofile->Close();


if (outTxtFile) fclose (outTxtFile);

// Now that everything is done, dump the list of lumiBlock that we processed in this job
if(!isMC){
	// FIXME: when lumi certificate is ready for rereco data, check that these work
	goodLumiFilter.FindLumiInFiles(urls);
	goodLumiFilter.DumpToJson(((outUrl.ReplaceAll(".root",""))+".json").Data());
	}

}

