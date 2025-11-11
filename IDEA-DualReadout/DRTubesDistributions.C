// A simple macro to plot the energy deposited, the S counts and the C counts and the light yields
// in the crystals from the IDEA_o2 k4geo simulation
// usage
// root -b
// .x DRTubesDistributions(path_to_file)

#include <iostream>
#include <string>

#include "podio/Reader.h"

#include "edm4hep/SimCalorimeterHitCollection.h"

int DRTubesDistributions(std::string input_file) {

  auto reader = podio::makeReader(input_file);

  TFile OutputFile = TFile("DRTubesDistributions.root", "RECREATE");
  TH1F DRBarrelS_th1 = TH1F("BarrelScintillationsignal", ";Scintillation [p.e.]; Number of events", 1000, 0.0, 4000.0);
  TH1F DRBarrelC_th1 = TH1F("BarrelCerenkovsignal", ";Cerenkov [p.e.]; Number of events", 1000, 0.0, 4000.0);
  TH1F DREndcapS_th1 = TH1F("EndcapScintillationsignal", ";Scintillation [p.e.]; Number of events", 1000, 0.0, 4000.0);
  TH1F DREndcapC_th1 = TH1F("EndcapCerenkovsignal", ";Cerenkov [p.e.]; Number of events", 1000, 0.0, 4000.0);

  // Fill the TH1 with the cell energy sum
  for (size_t i = 0; i < reader.getEvents(); ++i) {
    auto event = reader.readNextEvent();
    auto& BarrelS_hits = event.get<edm4hep::SimCalorimeterHitCollection>("DRBTScin");
    auto& BarrelC_hits = event.get<edm4hep::SimCalorimeterHitCollection>("DRBTCher");

    auto& EndcapLeftS_hits = event.get<edm4hep::SimCalorimeterHitCollection>("DRETScinLeft");
    auto& EndcapLeftC_hits = event.get<edm4hep::SimCalorimeterHitCollection>("DRETCherLeft");
    auto& EndcapRightS_hits = event.get<edm4hep::SimCalorimeterHitCollection>("DRETScinRight");
    auto& EndcapRightC_hits = event.get<edm4hep::SimCalorimeterHitCollection>("DRETCherRight");

    float TotalDRBTS = 0.;
    float TotalDRBTC = 0.;
    float TotalDRETS = 0.;
    float TotalDRETC = 0.;

    for (const auto& BarrelS_hit : BarrelS_hits) {
      TotalDRBTS += BarrelS_hit.getEnergy();
    }
    for (const auto& BarrelC_hit : BarrelC_hits) {
      TotalDRBTC += BarrelC_hit.getEnergy();
    }
    for (const auto& EndcapS_hit : EndcapLeftS_hits) {
      TotalDRETS += EndcapS_hit.getEnergy();
    }
    for (const auto& EndcapS_hit : EndcapRightS_hits) {
      TotalDRETS += EndcapS_hit.getEnergy();
    }
    for (const auto& EndcapC_hit : EndcapLeftC_hits) {
      TotalDRETC += EndcapC_hit.getEnergy();
    }
    for (const auto& EndcapC_hit : EndcapRightC_hits) {
      TotalDRETC += EndcapC_hit.getEnergy();
    }

    DRBarrelS_th1.Fill(TotalDRBTS);
    DRBarrelC_th1.Fill(TotalDRBTC);

    DREndcapS_th1.Fill(TotalDRETS);
    DREndcapC_th1.Fill(TotalDRETC);

    std::cout << "Barrel total S [p.e.] " << TotalDRBTS << " Barrel total C [p.e.] " << TotalDRBTC
              << " Endcap total S [p.e.] " << TotalDRETS << " Endcap total C [p.e.] " << TotalDRETC << std::endl;
  }

  DRBarrelS_th1.Write();
  DRBarrelC_th1.Write();

  DREndcapS_th1.Write();
  DREndcapC_th1.Write();

  OutputFile.Close();
  return 0;
}
