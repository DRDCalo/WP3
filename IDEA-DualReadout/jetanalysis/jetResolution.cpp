#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <numeric>
#include <cmath>
#include <sstream>

// PODIO includes
#include "podio/Reader.h"

// EDM4HEP includes
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/MutableCalorimeterHit.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/MCParticle.h"

// DD4HEP includes
#include "DDSegmentation/BitFieldCoder.h"

// FASTJET includes
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

// ROOT includes
#include "TH1F.h"
#include "TSystem.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TF1.h" 
#include "TStyle.h" 
#include "TFitResult.h"
#include "TPad.h"      // For ratio plots
#include "TLine.h"     // For ratio plots
#include "TLegend.h"   // For combined plot
#include "TPaveText.h" // For combined plot
#include "TMultiGraph.h"
#include "TLorentzVector.h"

// Function to get tower id from fiber id
int GetTowerNumber(long long int cellID, bool isBarrel){
    dd4hep::DDSegmentation::BitFieldCoder bcbarrel("system:5,stave:10,tower:-8,air:6,col:-16,row:16,clad:1,core:1,cherenkov:1");
    dd4hep::DDSegmentation:: BitFieldCoder bcendcap("system:5,stave:10,tower:6,air:1,col:16,row:16,clad:1,core:1,cherenkov:1");

    int tower;
    int stave;
    if (isBarrel){
        tower = bcbarrel.get(cellID,"tower");
	stave = bcbarrel.get(cellID,"stave");
    }
    else{
        tower = bcendcap.get(cellID,"tower");
	stave = bcendcap.get(cellID,"stave");
    }
    //std::cout<<" cellID "<<cellID<<" tower "<<tower<<" stave "<<stave<<std::endl;
    auto towerid = stave*80+tower;
    return towerid;
}

// Function to clusterize fiber hits into tower hits
edm4hep::CalorimeterHitCollection* aggregateToTowers(
    const edm4hep::SimCalorimeterHitCollection* inputHits,
    bool isBarrel
) {
    // The map uses the aggregated tower ID as key and stores the aggregated hit.
    // We use std::unordered_map for O(1) average lookup time.
    std::unordered_map<uint64_t, edm4hep::MutableCalorimeterHit> towerHitsMap;

    if (!inputHits) {
        std::cerr << "ERROR: Null inputHits pointer." << std::endl;
        return new edm4hep::CalorimeterHitCollection();
    }

    // 1. Iterate over every SimCaloHit (fiber) in the input collection
    for (const auto& simHit : *inputHits) {
        
        uint64_t fibreID = simHit.getCellID();
        
        // a) Get the aggregated Tower CellID
        uint64_t towerID = GetTowerNumber(fibreID, isBarrel);
        
        // b) Extract energy and position of the current fiber hit
        double hitEnergy = simHit.getEnergy();
        const auto& hitPos = simHit.getPosition();
        
        // c) Aggregation: check if the tower already exists in the map
        auto it = towerHitsMap.find(towerID);

        if (it == towerHitsMap.end()) {
            // FIRST HIT of the Tower: create a new CalorimeterHit and insert it
            edm4hep::MutableCalorimeterHit newHit;
            newHit.setCellID(towerID);
            newHit.setEnergy(hitEnergy);
            
            // Initial position is just the position of the first fiber
            newHit.setPosition(hitPos); 
            
            // Insert the hit into the map
            towerHitsMap.emplace(towerID, newHit);
        } else {
            // SUBSEQUENT HIT: accumulate energy and update position (weighted average)
            auto& towerHit = it->second;

            // Get current accumulated values of the tower
            double currentTowerEnergy = towerHit.getEnergy();
            const auto& currentTowerPos = towerHit.getPosition();

            // Calculate new total energy
            double newTotalEnergy = currentTowerEnergy + hitEnergy;

            // Calculate the new weighted position (Barycenter)
            // Formula: (P_old * E_old + P_new * E_new) / E_total
            if (newTotalEnergy > 0) {
                float weightedX = (currentTowerPos.x * currentTowerEnergy + hitPos.x * hitEnergy) / newTotalEnergy;
                float weightedY = (currentTowerPos.y * currentTowerEnergy + hitPos.y * hitEnergy) / newTotalEnergy;
                float weightedZ = (currentTowerPos.z * currentTowerEnergy + hitPos.z * hitEnergy) / newTotalEnergy;

                // Update the tower hit with the new position and energy
                towerHit.setPosition({weightedX, weightedY, weightedZ});
                towerHit.setEnergy(newTotalEnergy);
            } else {
                // Edge case: if energy is zero or negative (unlikely in sim), just add energy
                 towerHit.setEnergy(newTotalEnergy);
            }
        }
    } // End loop over fibers

    // 2. Create the output collection
    // NOTE: The caller is responsible for deallocating this pointer!
    edm4hep::CalorimeterHitCollection* outputCollection = new edm4hep::CalorimeterHitCollection();
    
    // 3. Fill the collection with the aggregated values from the map
    for (auto& pair : towerHitsMap) {
        // Use std::move to optimize transfer from map to vector
        outputCollection->push_back(std::move(pair.second));
    }
    
    return outputCollection;
};

// --- RESULTS STRUCTURE ---
// This struct holds the statistics extracted from the Gaussian fit of each histogram.
struct AnalysisResults {
    double meanMCTruth;
    double sigmaMCTruth;
    double semMCTruth;
    double sigmaMCTruthErr;
    
    double meanResidual;
    double sigmaResidual;
    double semResidual; 
    double sigmaResidualErr;
    
    double meanCombined;
    double sigmaCombined;
    double semCombined;
    double sigmaCombinedErr;
    
    // Pointers to the histograms
    TH1F* h_TotalDR; 
    TH1F* h_TotalMCTruth;
    TH1F* h_TotalResidual;
};

struct Angles {
    double theta;
    double phi;
};

Angles getThetaPhi(const edm4hep::CalorimeterHit& hit) {

    auto pos = hit.getPosition();
    double x = pos[0];
    double y = pos[1];
    double z = pos[2];

    double r_xy = std::sqrt(x * x + y * y);
    double phi = std::atan2(y, x);
    double theta = std::atan2(r_xy, z);
    return {theta, phi};
}

Angles getThetaPhi(const edm4hep::SimCalorimeterHit& hit) {

    auto pos = hit.getPosition();
    double x = pos[0];
    double y = pos[1];
    double z = pos[2];

    double r_xy = std::sqrt(x * x + y * y);
    double phi = std::atan2(y, x);
    double theta = std::atan2(r_xy, z);
    return {theta, phi};
}

double convertThetaToEta(double theta) {
    // eta = -ln(tan(theta/2))
    return -1.0 * std::log(std::tan(theta / 2.0));
}

fastjet::PseudoJet mergejet(fastjet::PseudoJet jet_scin, fastjet::PseudoJet jet_cher, double chi) {

  double jetPx = (jet_scin.px()-chi*jet_cher.px())/(1-chi);
  double jetPy = (jet_scin.py()-chi*jet_cher.py())/(1-chi);
  double jetPz = (jet_scin.pz()-chi*jet_cher.pz())/(1-chi);
  double jetE = (jet_scin.e()-chi*jet_cher.e())/(1.-chi);
  return fastjet::PseudoJet(jetPx, jetPy, jetPz, jetE);
}

fastjet::PseudoJet matchjet(fastjet::PseudoJet jet_in, std::vector<fastjet::PseudoJet> testvec) {

  int imin=-1;
  double deltarmin=99999.;
  for(uint i=0; i<testvec.size(); i++){
    double deltar=jet_in.delta_R(testvec.at(i));
    if(deltar<deltarmin) {
      deltarmin=deltar;
      imin=i;
    }
  }
  if(imin != -1) return testvec.at(imin);
  else
  return fastjet::PseudoJet(0., 0., 0., 0.);
}

/**
 * @brief Analyzes a single input file, creates histograms in GeV, and fits them.
 * @param input_file Path to the PODIO file.
 * @param energy The beam energy (GeV), used for histogram naming.
 * @param Ks Scintillation calibration constant (p.e./GeV).
 * @param Kc Cerenkov calibration constant (p.e./GeV).
 * @return An AnalysisResults struct with stats from Gaussian fits and histogram pointers.
 */
AnalysisResults analyzeFile(const std::string& input_file, double energy, double Ks, double Kc, double Chi, double Ks_crystal, double Kc_crystal, double Chi_crystal) {
    using namespace dd4hep;
 
    std::cout << "  Opening file: " << input_file << std::endl;
    auto reader = podio::makeReader(input_file);

    // --- Create unique names and titles for histograms (all in GeV) ---
    std::stringstream ss_s_name, ss_c_name, ss_e_name, ss_s_title, ss_c_title, ss_e_title;
    ss_s_name << "h_TotalDR_E_" << energy;
    ss_c_name << "h_TotalMCTruth_E_" << energy;
    ss_e_name << "h_Residual_E_" << energy;
    
    // Update titles to reflect GeV units
    ss_s_title << "Energy DR (E=" << energy << " GeV);DR jet energy [GeV]; Events";
    ss_c_title << "Energy MC truth (E=" << energy << " GeV);MCTruth jet energy [GeV]; Events";
    ss_e_title << "Residual Energy (E=" << energy << " GeV);(Ej-Emc)/Emc; Events";

    // Create histograms with 'new' and SetDirectory(nullptr) to make them persist.
    // Update ranges to be appropriate for GeV.
    double range_min = 0.0;
    double range_max = energy * 1.5;
    TH1F* h_TotalDR = new TH1F(ss_s_name.str().c_str(), ss_s_title.str().c_str(), 200, range_min, range_max); 
    TH1F* h_TotalMCTruth = new TH1F(ss_c_name.str().c_str(), ss_c_title.str().c_str(), 200, range_min, range_max);
    TH1F* h_TotalResidual = new TH1F(ss_e_name.str().c_str(), ss_e_title.str().c_str(), 200, -2, 2);

    h_TotalDR->SetDirectory(nullptr);
    h_TotalMCTruth->SetDirectory(nullptr);
    h_TotalResidual->SetDirectory(nullptr);

    unsigned int nEvents = reader.getEvents();
    std::cout << "  Number of events in file: " << nEvents << std::endl;

    // --- Event Loop ---
    for (size_t i = 0; i < nEvents; ++i) {

        auto event = reader.readNextEvent();

        const auto& mcParticles = event.get<edm4hep::MCParticleCollection>("MCParticles");

	// Get hadronic calorimeter fiber hit collections
        auto& BarrelS_hits = event.get<edm4hep::SimCalorimeterHitCollection>("DRBTScin");
        auto& BarrelC_hits = event.get<edm4hep::SimCalorimeterHitCollection>("DRBTCher");

        auto& EndcapLeftS_hits = event.get<edm4hep::SimCalorimeterHitCollection>("DRETScinLeft");
        auto& EndcapLeftC_hits = event.get<edm4hep::SimCalorimeterHitCollection>("DRETCherLeft");
        auto& EndcapRightS_hits = event.get<edm4hep::SimCalorimeterHitCollection>("DRETScinRight");
        auto& EndcapRightC_hits = event.get<edm4hep::SimCalorimeterHitCollection>("DRETCherRight");
        
	// Aggegate hadronic calorimeter fiber collections to tower collections
        auto STowerBarrelCollection = aggregateToTowers(&BarrelS_hits, true);
        auto CTowerBarrelCollection = aggregateToTowers(&BarrelC_hits, true);
        auto STowerEndcapLCollection = aggregateToTowers(&EndcapLeftS_hits, false);
        auto CTowerEndcapLCollection = aggregateToTowers(&EndcapLeftC_hits, false);
        auto STowerEndcapRCollection = aggregateToTowers(&EndcapRightS_hits, false);
        auto CTowerEndcapRCollection = aggregateToTowers(&EndcapRightC_hits, false);

	// Get crystal hit collections
        auto& crystals_Scounts = event.get<edm4hep::SimCalorimeterHitCollection>("SCEPCal_MainScounts");
        auto& crystals_Ccounts = event.get<edm4hep::SimCalorimeterHitCollection>("SCEPCal_MainCcounts");

	// Collection of inputs to fastjet
	std::vector<fastjet::PseudoJet> inputparticles_scin;
	std::vector<fastjet::PseudoJet> inputparticles_cher;
	std::vector<fastjet::PseudoJet> inputparticles_scin_crystal;
	std::vector<fastjet::PseudoJet> inputparticles_cher_crystal;
	std::vector<fastjet::PseudoJet> inputparticles_mctruth;
	std::vector<fastjet::PseudoJet> jet_scin;
	std::vector<fastjet::PseudoJet> jet_cher;
	std::vector<fastjet::PseudoJet> jet_cher_aligned;
	std::vector<fastjet::PseudoJet> jet_scin_crystal;
	std::vector<fastjet::PseudoJet> jet_cher_crystal;
	std::vector<fastjet::PseudoJet> jet_cher_crystal_aligned;
	std::vector<fastjet::PseudoJet> jet_dr;
	std::vector<fastjet::PseudoJet> jet_dr_crystal;
	std::vector<fastjet::PseudoJet> jet_dr_crystal_aligned;
	std::vector<fastjet::PseudoJet> jet_final;
	std::vector<fastjet::PseudoJet> jet_mctruth;
	std::vector<fastjet::PseudoJet> jet_mctruth_aligned;

	// Create input objects for fastjet
        for (const auto& anhit : *STowerBarrelCollection) { // scintillating hadronic barrel
	    double energy_calibrated = anhit.getEnergy() / Ks;
	    auto angles = getThetaPhi(anhit);
            double pt = energy_calibrated*sin(angles.theta);
            TLorentzVector tower;
            tower.SetPtEtaPhiM(pt, convertThetaToEta(angles.theta), angles.phi, 0.);
	    inputparticles_scin.push_back( fastjet::PseudoJet(tower.Px(), tower.Py(), tower.Pz(), tower.E()) );;
	}
        for (const auto& anhit : *STowerEndcapLCollection) { // scintillating hadronic endcap left
	    double energy_calibrated = anhit.getEnergy() / Ks;
	    auto angles = getThetaPhi(anhit);
            double pt = energy_calibrated*sin(angles.theta);
            TLorentzVector tower;
            tower.SetPtEtaPhiM(pt, convertThetaToEta(angles.theta), angles.phi, 0.);
	    inputparticles_scin.push_back( fastjet::PseudoJet(tower.Px(), tower.Py(), tower.Pz(), tower.E()) );;
	}
        for (const auto& anhit : *STowerEndcapRCollection) { // scintillating hadronic endcap right
	    double energy_calibrated = anhit.getEnergy() / Ks;
	    auto angles = getThetaPhi(anhit);
            double pt = energy_calibrated*sin(angles.theta);
            TLorentzVector tower;
            tower.SetPtEtaPhiM(pt, convertThetaToEta(angles.theta), angles.phi, 0.);
	    inputparticles_scin.push_back( fastjet::PseudoJet(tower.Px(), tower.Py(), tower.Pz(), tower.E()) );;
	}

	fastjet::JetDefinition jet_defs(fastjet::ee_genkt_algorithm, 2.*M_PI, 1.);
        fastjet::ClusterSequence clust_seq_scin(inputparticles_scin, jet_defs); 
        jet_scin.push_back(clust_seq_scin.exclusive_jets(int(2))[0]);
        jet_scin.push_back(clust_seq_scin.exclusive_jets(int(2))[1]);

        for (const auto& anhit : crystals_Scounts) { // scintillating crystals
	    double energy_calibrated = anhit.getEnergy() / Ks_crystal;
	    auto angles = getThetaPhi(anhit);
            double pt = energy_calibrated*sin(angles.theta);
	    if(energy_calibrated > 0.0001) { // GeV
              TLorentzVector tower;
              tower.SetPtEtaPhiM(pt, convertThetaToEta(angles.theta), angles.phi, 0.);
	      inputparticles_scin_crystal.push_back( fastjet::PseudoJet(tower.Px(), tower.Py(), tower.Pz(), tower.E()) );;
	    }
	}

        fastjet::ClusterSequence clust_seq_scin_crystal(inputparticles_scin_crystal, jet_defs); 
        jet_scin_crystal.push_back(clust_seq_scin_crystal.exclusive_jets(int(2))[0]);
        jet_scin_crystal.push_back(clust_seq_scin_crystal.exclusive_jets(int(2))[1]);

        for (const auto& anhit : *CTowerBarrelCollection) { // Cerenkov hadronic barrel
	    double energy_calibrated = anhit.getEnergy() / Kc;
	    auto angles = getThetaPhi(anhit);
            double pt = energy_calibrated*sin(angles.theta);
            TLorentzVector tower;
            tower.SetPtEtaPhiM(pt, convertThetaToEta(angles.theta), angles.phi, 0.);
	    inputparticles_cher.push_back( fastjet::PseudoJet(tower.Px(), tower.Py(), tower.Pz(), tower.E()) );;
	}
        for (const auto& anhit : *CTowerEndcapLCollection) { // Cerenkov hadronic endcap left
	    double energy_calibrated = anhit.getEnergy() / Kc;
	    auto angles = getThetaPhi(anhit);
            double pt = energy_calibrated*sin(angles.theta);
            TLorentzVector tower;
            tower.SetPtEtaPhiM(pt, convertThetaToEta(angles.theta), angles.phi, 0.);
	    inputparticles_cher.push_back( fastjet::PseudoJet(tower.Px(), tower.Py(), tower.Pz(), tower.E()) );;
	}
        for (const auto& anhit : *CTowerEndcapRCollection) { // Cerenkov hadronic endcap right
	    double energy_calibrated = anhit.getEnergy() / Kc;
	    auto angles = getThetaPhi(anhit);
            double pt = energy_calibrated*sin(angles.theta);
            TLorentzVector tower;
            tower.SetPtEtaPhiM(pt, convertThetaToEta(angles.theta), angles.phi, 0.);
	    inputparticles_cher.push_back( fastjet::PseudoJet(tower.Px(), tower.Py(), tower.Pz(), tower.E()) );;
	}

        fastjet::ClusterSequence clust_seq_cher(inputparticles_cher, jet_defs); 
        jet_cher.push_back(clust_seq_cher.exclusive_jets(int(2))[0]);
        jet_cher.push_back(clust_seq_cher.exclusive_jets(int(2))[1]);
        
	for (const auto& anhit : crystals_Ccounts) { // Cerenkov crystals
	    double energy_calibrated = anhit.getEnergy() / Kc_crystal;
	    auto angles = getThetaPhi(anhit);
            double pt = energy_calibrated*sin(angles.theta);
	    if(energy_calibrated > 0.0001) { // GeV
              TLorentzVector tower;
              tower.SetPtEtaPhiM(pt, convertThetaToEta(angles.theta), angles.phi, 0.);
	      inputparticles_cher_crystal.push_back( fastjet::PseudoJet(tower.Px(), tower.Py(), tower.Pz(), tower.E()) );;
	    }
	}

        fastjet::ClusterSequence clust_seq_cher_crystal(inputparticles_cher_crystal, jet_defs); 
        jet_cher_crystal.push_back(clust_seq_cher_crystal.exclusive_jets(int(2))[0]);
        jet_cher_crystal.push_back(clust_seq_cher_crystal.exclusive_jets(int(2))[1]);

	/*for(auto SPseudoJet : jet_scin){
	    std::cout<<"pseudojet scin eta "<<SPseudoJet.eta()<<" energy "<<SPseudoJet.E()<<std::endl;
	};
	for(auto CPseudoJet : jet_cher){
	    std::cout<<"pseudojet cher eta "<<CPseudoJet.eta()<<" energy "<<CPseudoJet.E()<<std::endl;
	};*/
        
	// Align cher jets to scin jets
	for(uint jn=0; jn<jet_scin.size(); jn++) {
            jet_cher_aligned.push_back(matchjet(jet_scin[jn], jet_cher)); 
        }
	// Align crystal cher jets to scin jets
	for(uint jn=0; jn<jet_scin_crystal.size(); jn++) {
            jet_cher_crystal_aligned.push_back(matchjet(jet_scin_crystal[jn], jet_cher_crystal)); 
        }

	/*std::cout<<"After alignment"<<std::endl;
	for(auto SPseudoJet : jet_scin){
	    std::cout<<"pseudojet scin eta "<<SPseudoJet.eta()<<" energy "<<SPseudoJet.E()<<std::endl;
	};
	for(auto CPseudoJet : jet_cher_aligned){
	    std::cout<<"pseudojet cher eta "<<CPseudoJet.eta()<<" energy "<<CPseudoJet.E()<<std::endl;
	};*/
       
	// Apply dual-readout correction to jets
        for(uint jn=0; jn<jet_scin.size(); jn++) {
            jet_dr.push_back(mergejet(jet_scin[jn],jet_cher_aligned[jn], Chi));
        }

	// Apply dual-readout correction to crystal jets
        for(uint jn=0; jn<jet_scin_crystal.size(); jn++) {
            jet_dr_crystal.push_back(mergejet(jet_scin_crystal[jn],jet_cher_crystal_aligned[jn], Chi_crystal));
        }

	/*for(auto DRPseudoJet : jet_dr){
	    std::cout<<"pseudojet dr eta "<<DRPseudoJet.eta()<<" energy "<<DRPseudoJet.E()<<std::endl;
	};
	for(auto DRPseudoJet : jet_dr_crystal){
	    std::cout<<"pseudojet crystal dr eta "<<DRPseudoJet.eta()<<" energy "<<DRPseudoJet.E()<<std::endl;
	};*/

	// Align crystal jets to fiber jets
	for(uint jn=0; jn<jet_dr.size(); jn++) {
            jet_dr_crystal_aligned.push_back(matchjet(jet_dr[jn], jet_dr_crystal)); 
        }

	// Sum crystal and fiber aligned jets
	for(uint jn=0; jn<jet_dr.size(); jn++) {
	    jet_final.push_back(jet_dr[jn]+jet_dr_crystal_aligned[jn]);
	}

	/*for(auto PseudoJet : jet_final){
	    std::cout<<"pseudojet final eta "<<PseudoJet.eta()<<" energy "<<PseudoJet.E()<<std::endl;
	};*/

        for (const auto& p : mcParticles) { // MC truth particles
	    if(p.getGeneratorStatus() != 1) continue;
	    auto pdgid = p.getPDG();
	    if(std::abs(pdgid) == 13 || std::abs(pdgid) == 12 || std::abs(pdgid) == 14 || std::abs(pdgid)== 16) continue;
	    //auto daughters = p.getDaughters();
	    //if(daughters.size() > 0) continue;
	    //std::cout<<" pdg "<<p.getPDG()<<" daughters "<<daughters.size()<<" status "<<p.getGeneratorStatus()<<std::endl;
	    double E = p.getEnergy();
	    auto Momentum = p.getMomentum();
            float px = Momentum.x;
            float py = Momentum.y;
            float pz = Momentum.z;
	    inputparticles_mctruth.push_back(fastjet::PseudoJet(px, py, pz, E));
	}

        fastjet::ClusterSequence clust_seq_mctruth(inputparticles_mctruth, jet_defs); 
        jet_mctruth.push_back(clust_seq_mctruth.exclusive_jets(int(2))[0]);
        jet_mctruth.push_back(clust_seq_mctruth.exclusive_jets(int(2))[1]);

	/*for(auto MCPseudoJet : jet_mctruth){
	    std::cout<<"mcpseudojet eta "<<MCPseudoJet.eta()<<" energy "<<MCPseudoJet.E()<<std::endl;
	};*/

	// Align final jets to mc truth jets
	for(uint jn=0; jn<jet_final.size(); jn++) {
            jet_mctruth_aligned.push_back(matchjet(jet_final[jn], jet_mctruth));
        }

	for(std::size_t t=0; t<2; t++){
	    std::cout<<"reco jet "<<t<<" eta "<<jet_final[t].eta()<<" phi "<<jet_final[t].phi()<<" energy "<<jet_final[t].E()<<" mc truth jet "<<t<<" eta "<<jet_mctruth_aligned[t].eta()<<" phi "<<jet_mctruth_aligned[t].phi()<<" energy "<<jet_mctruth_aligned[t].E()<<std::endl;
	}

	// Fill histograms with calibrated GeV values
	if( jet_mctruth_aligned[0].eta()>-1.74 && jet_mctruth_aligned[0].eta()<1.74){
            h_TotalDR->Fill(jet_final[0].E());
            h_TotalMCTruth->Fill(jet_mctruth_aligned[0].E());
            h_TotalResidual->Fill((jet_final[0].E()-jet_mctruth_aligned[0].E())/jet_mctruth_aligned[0].E());
	}
	std::cout<<"end of event"<<std::endl;
    } // --- End Event Loop ---
    
    // --- Calculate Statistics from Gaussian Fit ---
    AnalysisResults result;

    // Helper function to fit a histogram and extract parameters
    auto fitHistogram = [&](TH1F* h, double energyfit) -> std::tuple<double, double, double, double> {
        if (h->GetEntries() == 0) {
            std::cerr << "  WARNING: Histogram " << h->GetName() << " is empty. Skipping fit." << std::endl;
            return {0, 0, 0, 0};
        }
        
        // Initial parameter estimates from histogram
        double initialMean = h->GetMean();
        double initialSigma = h->GetRMS();
        double fitMin = initialMean - 1.0 * initialSigma;
        double fitMax = initialMean + 1.0 * initialSigma;

	if(energyfit==20){
            fitMin = -0.22;
            fitMax = 0.04;
	}

        TF1* gausFit = new TF1("gausFit", "gaus", fitMin, fitMax);
        gausFit->SetParameters(h->GetMaximum(), initialMean, initialSigma);
        
        // Perform the fit ("R" = use range, "Q" = quiet, "S" = get TFitResultPtr)
        TFitResultPtr fitResult = h->Fit(gausFit, "RQS"); 
        
        if (!fitResult.Get() || !fitResult->IsValid() || (int)fitResult != 0) {
             std::cerr << "  WARNING: Gaussian fit failed for histogram " << h->GetName() << ". Using GetMean/GetRMS." << std::endl;
             delete gausFit;
             double mean = h->GetMean();
             double sigma = h->GetRMS();
             double sem = (h->GetEntries() > 0) ? sigma / std::sqrt(h->GetEntries()) : 0;
             double serr = (h->GetEntries() > 1) ? sigma / std::sqrt(2.0 * (h->GetEntries() - 1)) : 0; // Approx error on sigma
             return {mean, sigma, sem, serr};
        }

        double mean = fitResult->Parameter(1);
        double meanErr = fitResult->ParError(1);
        double sigma = fitResult->Parameter(2);
        double sigmaErr = fitResult->ParError(2);
        
        delete gausFit;
        return {mean, sigma, meanErr, sigmaErr};
    };

    auto [meanMCTruth, sigmaMCTruth, semMCTruth, sigmaMCTruthErr] = fitHistogram(h_TotalMCTruth, energy);
    result.meanMCTruth = meanMCTruth;
    result.sigmaMCTruth = sigmaMCTruth;
    result.semMCTruth = semMCTruth;
    result.sigmaMCTruthErr = sigmaMCTruthErr;

    auto [meanC, sigmaC, semC, sigmaCErr] = fitHistogram(h_TotalDR, energy);
    result.meanCombined = meanC;
    result.sigmaCombined = sigmaC;
    result.semCombined = semC;
    result.sigmaCombinedErr = sigmaCErr;

    auto [meanE, sigmaE, semE, sigmaEErr] = fitHistogram(h_TotalResidual, energy);
    result.meanResidual = meanE;
    result.sigmaResidual = sigmaE;
    result.semResidual = semE;
    result.sigmaResidualErr = sigmaEErr;

    // Store histogram pointers
    result.h_TotalMCTruth = h_TotalMCTruth;
    result.h_TotalDR = h_TotalDR;
    result.h_TotalResidual = h_TotalResidual;

    return result;
}

/**
 * @brief Main function to run the analysis and produce linearity and resolution graphs.
 */
int main() {
    
    // --- CALIBRATION CONSTANTS (p.e./GeV) ---
    const double CALIBRATION_CONSTANT_S = 206.284; // S p.e./GeV
    const double CALIBRATION_CONSTANT_C = 68.069; // C p.e./GeV
    const double CALIBRATION_CONSTANT_Chi = 0.31; // hadronic scale calibration constant
    
    // --- CALIBRATION CONSTANTS CRYSTALS (p.e./GeV) ---
    const double CALIBRATION_CRYSTAL_S = 1960.92; // S p.e./GeV
    const double CALIBRATION_CRYSTAL_C = 97.531; // C p.e./GeV
    const double CALIBRATION_CRYSTAL_Chi = 0.41; // hadronic scale calibration constant for crystals

    // --- FILE MAP (Energy in GeV) ---
    std::map<double, std::string> filesToAnalyze;
    filesToAnalyze[20.0] = "../IDEA_o2_v01_jet20gev.root";
    filesToAnalyze[40.0] = "../IDEA_o2_v01_jet40gev.root";
    filesToAnalyze[60.0] = "../IDEA_o2_v01_jet60gev.root";
    // ... add more energies and files here

    std::vector<TH1F*> allHistograms;
    std::vector<TObject*> objectsToClean; // Store all created ROOT objects

    /*TGraphErrors* grE_Linearity = new TGraphErrors();
    grE_Linearity->SetName("g_Linearity_Combined");
    grE_Linearity->SetTitle("Combined Signal Linearity;True Energy [GeV];Dual-readout energy [GeV]");
    grE_Linearity->SetMarkerStyle(20); grE_Linearity->SetMarkerColor(kGreen+2); grE_Linearity->SetLineColor(kGreen+2);
    objectsToClean.push_back(grE_Linearity);*/

    // --- 2. Resolution Graphs (Sigma/Mean vs. Energy) ---
    TGraphErrors* grE_Resolution = new TGraphErrors();
    grE_Resolution->SetName("g_Resolution_Combined");
    grE_Resolution->SetTitle("Dual-readout Jet Energy Resolution (crystal+fibers);E_{CM}/2;Resolution #sigma((E_{j}-E_{mc})/E_{mc})");
    grE_Resolution->SetMarkerStyle(20); grE_Resolution->SetMarkerColor(kGreen+2); grE_Resolution->SetLineColor(kGreen+2);
    objectsToClean.push_back(grE_Resolution);

    int pointIndex = 0; 
    double minEnergy = 1e9, maxEnergy = -1;

    std::cout << "Starting energy scan analysis with calibration..." << std::endl;
    std::cout << "Calibration S: " << CALIBRATION_CONSTANT_S << " p.e./GeV" << std::endl;
    std::cout << "Calibration C: " << CALIBRATION_CONSTANT_C << " p.e./GeV" << std::endl;
    std::cout << "Calibration crystal S: " << CALIBRATION_CRYSTAL_S << " p.e./GeV" << std::endl;
    std::cout << "Calibration crystal C: " << CALIBRATION_CRYSTAL_C << " p.e./GeV" << std::endl;
    std::cout << "Hadronic Chi factor: " << CALIBRATION_CONSTANT_Chi << std::endl;
    std::cout << "Crystal Chi factor: " << CALIBRATION_CRYSTAL_Chi << std::endl;

    for (const auto& pair : filesToAnalyze) {
        double energy = pair.first;
        std::string filename = pair.second;

        if (energy/2. < minEnergy) minEnergy = energy/2.-1.;
        if (energy/2. > maxEnergy) maxEnergy = energy/2.+1.;

        if (gSystem->AccessPathName(filename.c_str())) {
            std::cerr << "ERROR: Cannot find file: " << filename << std::endl;
            continue; 
        }

        std::cout << "\nAnalyzing for energy = " << energy << " GeV (File: " << filename << ")" << std::endl;
        
        AnalysisResults results = analyzeFile(filename, energy, CALIBRATION_CONSTANT_S, CALIBRATION_CONSTANT_C, CALIBRATION_CONSTANT_Chi, CALIBRATION_CRYSTAL_S, CALIBRATION_CRYSTAL_C, CALIBRATION_CRYSTAL_Chi);

        if (results.meanMCTruth > 0) {
            
            double resE = results.sigmaResidual;
	    double err_resE = results.sigmaResidualErr;

            grE_Resolution->SetPoint(pointIndex, energy/2., resE);
            grE_Resolution->SetPointError(pointIndex, 0.0, err_resE);

            pointIndex++; 

            // Store histogram pointers for saving
            if (results.h_TotalMCTruth) allHistograms.push_back(results.h_TotalMCTruth);
            if (results.h_TotalDR) allHistograms.push_back(results.h_TotalDR);
            if (results.h_TotalResidual) allHistograms.push_back(results.h_TotalResidual);

        } else {
            std::cerr << "  No valid results from analysis for energy = " << energy << " GeV." << std::endl;
            delete results.h_TotalMCTruth;
            delete results.h_TotalDR;
            delete results.h_TotalResidual;
        }
    }

    // --- Create and save TGraphErrors, Fits, and Histograms ---
    if (pointIndex > 0) {
        std::cout << "\nSaving Calibrated Graphs, Fits, and Histograms to 'jetEnergyScanCalibratedResults.root'..." << std::endl;
        
        gStyle->SetOptFit(0); // Turn off global fit stats box
        TFile *outFile = new TFile("jetEnergyScanCalibratedResults.root", "RECREATE");

        // --- Resolution ---
        outFile->mkdir("Resolution");
        outFile->cd("Resolution");
        
        TF1* fit_res_E = new TF1("fit_res_E", "sqrt( ([0]/sqrt(x))^2 + [1]^2 )", minEnergy, maxEnergy);
        fit_res_E->SetParNames("Stochastic", "Constant");
        fit_res_E->SetParameters(0.1, 0.01); 
        fit_res_E->SetLineColor(kGreen+2);
        grE_Resolution->Fit(fit_res_E, "RQ");
        objectsToClean.push_back(fit_res_E);

        // E Combined Resolution Canvas
        TCanvas *cE_res = new TCanvas("c_Resolution_Combined", "Energy Resolution", 800, 600);
        cE_res->SetLeftMargin(0.15); // Fix Y label cutoff
        cE_res->SetGrid();
        grE_Resolution->Draw("AP");
        fit_res_E->Draw("SAME");
        // PaveText for S Fit Results
        TPaveText *ptDRE = new TPaveText(0.18, 0.75, 0.53, 0.88, "NDC");
        ptDRE->SetBorderSize(1);
        ptDRE->SetFillColor(0);
        ptDRE->SetTextColor(kGreen+2);
        ptDRE->SetTextAlign(12); // Left-aligned
        ptDRE->AddText(Form("DR Fit: a = %.3f #pm %.4f", fit_res_E->GetParameter(0), fit_res_E->GetParError(0)));
        ptDRE->AddText(Form("DR Fit: c = %.3f #pm %.4f", fit_res_E->GetParameter(1), fit_res_E->GetParError(1)));
        ptDRE->Draw();
        cE_res->Write();
        grE_Resolution->Write();
        objectsToClean.push_back(ptDRE);
        objectsToClean.push_back(cE_res);

        // --- Write Histograms ---
        // All histograms are now in GeV
        outFile->mkdir("Histograms_GeV");
        outFile->cd("Histograms_GeV");
        for (auto h : allHistograms) { if (h) h->Write(); }

        outFile->Close();
        std::cout << "Analysis results saved to 'jetEnergyScanCalibratedResults.root'." << std::endl;
        
        // --- Clean up memory ---
        for (auto h : allHistograms) delete h;
        allHistograms.clear();
        
        // Delete all other ROOT objects created with 'new'
        // TCanvases are cleaned up by TFile::Close() if they were written
        // but graphs, fits, legends, etc. are not.
        for (auto obj : objectsToClean) delete obj;
        objectsToClean.clear();

    } else {
        std::cout << "\nNo valid points to draw." << std::endl;
        // Clean up empty graphs and histograms
        for (auto obj : objectsToClean) delete obj;
        objectsToClean.clear();
        for (auto h : allHistograms) delete h;
    }

    return 0;
}
