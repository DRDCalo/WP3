// A simple macro to plot the energy deposited, the S counts and the C counts and the light yields
// in the crystals from the IDEA_o2 k4geo simulation
// usage
// root -b
// .x SumCrystalsInformation(path_to_file)

#include <iostream>
#include <string>

#include "podio/Reader.h"

#include "edm4hep/SimCalorimeterHitCollection.h"

int SumCrystalsInformation(std::string input_file) {

  auto reader = podio::makeReader(input_file);

  TFile OutputFile = TFile("SumCrystalsInformation.root", "RECREATE");
  TH1F EnergyDeposited_th1 = TH1F("Energy deposited", ";energy sum [GeV]; Number of events", 1000, 70, 12.0);
  TH1F Scounts_th1 = TH1F("Scintillation counts", ";Scintillation counts; Number of events", 1000, 0, 1000000.0);
  TH1F Ccounts_th1 = TH1F("Cerenkov counts", ";Cerenkov counts sum; Number of events", 1000, 0, 1000000.0);
  TH1F SLY_th1 = TH1F("S light yield", ";Scintillation light yield [p.e./GeV]; Number of events", 1000, 0, 100000.0);
  TH1F CLY_th1 = TH1F("C light yield", ";Cerenkov light yield [p.e./GeV]; Number of events", 1000, 0, 100000.0);

  // Fill the TH1 with the cell energy sum
  for (size_t i = 0; i < reader.getEvents(); ++i) {
    auto event = reader.readNextEvent();
    auto& crystals_edep = event.get<edm4hep::SimCalorimeterHitCollection>("SCEPCal_MainEdep");
    auto& crystals_Scounts = event.get<edm4hep::SimCalorimeterHitCollection>("SCEPCal_MainScounts");
    auto& crystals_Ccounts = event.get<edm4hep::SimCalorimeterHitCollection>("SCEPCal_MainCcounts");
    float total_energy = 0.;
    float Scounts = 0.;
    float Ccounts = 0.;
    for (const auto& crystal_edep : crystals_edep) {
      total_energy += crystal_edep.getEnergy();
    }
    for (const auto& S : crystals_Scounts) {
      Scounts += S.getEnergy();
    }
    for (const auto& C : crystals_Ccounts) {
      Ccounts += C.getEnergy();
    }
    EnergyDeposited_th1.Fill(total_energy);
    Scounts_th1.Fill(Scounts);
    Ccounts_th1.Fill(Ccounts);
    SLY_th1.Fill(Scounts / total_energy);
    CLY_th1.Fill(Ccounts / total_energy);
    std::cout << "Energy deposited in crystals [GeV]: " << total_energy << " S counts " << Scounts << " C counts "
              << Ccounts << std::endl;
  }
  EnergyDeposited_th1.Write();
  Scounts_th1.Write();
  Ccounts_th1.Write();
  SLY_th1.Write();
  CLY_th1.Write();

  OutputFile.Close();
  return 0;
}
