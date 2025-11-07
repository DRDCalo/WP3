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
  TH1F EnergyDeposited_th1 = TH1F("Energy deposited", ";energy sum [GeV]; Number of events", 1000, 7, 12.0);
  TH1F Scounts_th1 = TH1F("Scintillation counts", ";Scintillation counts; Number of events", 1000, 0, 40000.0);
  TH1F Ccounts_th1 = TH1F("Cerenkov counts", ";Cerenkov counts sum; Number of events", 100, 0, 2000.0);
  TH1F SLY_th1 = TH1F("S light yield", ";Scintillation light yield [p.e./GeV]; Number of events", 100, 0, 4000.0);
  TH1F CLY_th1 = TH1F("C light yield", ";Cerenkov light yield [p.e./GeV]; Number of events", 100, 0, 200.0);

  TH1F TimingLayerEdep_th1 =
      TH1F("Timing layer Energy deposited", ";energy sum [GeV]; Number of events", 100, 0.0, 0.5);
  TH1F TimingLayerScounts_th1 =
      TH1F("Timing layer Scintillation counts", ";Scintillation counts; Number of events", 500, 0, 200000.0);
  TH1F TimingLayerSLY_th1 =
      TH1F("Timing layer S light yield", ";Scintillation light yield [p.e./MeV]; Number of events", 1000, 0, 10000.0);

  // Fill the TH1 with the cell energy sum
  for (size_t i = 0; i < reader.getEvents(); ++i) {
    auto event = reader.readNextEvent();
    auto& crystals_edep = event.get<edm4hep::SimCalorimeterHitCollection>("SCEPCal_MainEdep");
    auto& crystals_Scounts = event.get<edm4hep::SimCalorimeterHitCollection>("SCEPCal_MainScounts");
    auto& crystals_Ccounts = event.get<edm4hep::SimCalorimeterHitCollection>("SCEPCal_MainCcounts");
    auto& TimingLayer_edeps = event.get<edm4hep::SimCalorimeterHitCollection>("SCEPCal_TimingEdep");
    auto& TimingLayer_Scounts = event.get<edm4hep::SimCalorimeterHitCollection>("SCEPCal_TimingScounts");

    float total_energy = 0.;
    float Scounts = 0.;
    float Ccounts = 0.;
    float TL_total_energy = 0.;
    float TL_Scounts = 0.;
    for (const auto& crystal_edep : crystals_edep) {
      total_energy += crystal_edep.getEnergy();
    }
    for (const auto& S : crystals_Scounts) {
      Scounts += S.getEnergy();
    }
    for (const auto& C : crystals_Ccounts) {
      Ccounts += C.getEnergy();
    }
    for (const auto& tl_edep : TimingLayer_edeps) {
      TL_total_energy += tl_edep.getEnergy();
    }
    for (const auto& tl_S : TimingLayer_Scounts) {
      TL_Scounts += tl_S.getEnergy();
    }

    EnergyDeposited_th1.Fill(total_energy);
    Scounts_th1.Fill(Scounts);
    Ccounts_th1.Fill(Ccounts);
    SLY_th1.Fill(Scounts / total_energy);
    CLY_th1.Fill(Ccounts / total_energy);

    TimingLayerEdep_th1.Fill(TL_total_energy);
    TimingLayerScounts_th1.Fill(TL_Scounts);
    TimingLayerSLY_th1.Fill(TL_Scounts / (TL_total_energy * 1000)); // LY expressedn in S counts / MeV

    std::cout << "Energy deposited in crystals [GeV]: " << total_energy << " S counts " << Scounts << " C counts "
              << Ccounts << " Energy deposited in timing layer [GeV] " << TL_total_energy << " S counts " << TL_Scounts
              << std::endl;
  }
  EnergyDeposited_th1.Write();
  Scounts_th1.Write();
  Ccounts_th1.Write();
  SLY_th1.Write();
  CLY_th1.Write();
  TimingLayerEdep_th1.Write();
  TimingLayerScounts_th1.Write();
  TimingLayerSLY_th1.Write();

  OutputFile.Close();
  return 0;
}
