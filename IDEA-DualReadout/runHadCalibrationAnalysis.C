// A simple macro to extract calibration constants for the tubes calorimeter at different theta values
#include <cmath> // For std::sqrt
#include <iostream>
#include <map>
#include <numeric> // For std::accumulate
#include <sstream> // For unique histogram names
#include <string>
#include <vector>

#include "podio/Reader.h"
//#include "podio/podio-config.h" // For PODIO_VERSION

#include "edm4hep/SimCalorimeterHitCollection.h"

#include "TCanvas.h"      // Added to save canvases
#include "TFile.h"        // Added to save graphs
#include "TGraphErrors.h" // Added for graphs
#include "TH1F.h"
#include "TSystem.h"

// Structure to hold the analysis results for a single file
struct AnalysisResults {
  double meanScint;     // Mean of the total S signal (Barrel + Endcap)
  double sigmaScint;    // Sigma (RMS) of the total S signal
  double semScint;      // Standard Error on the Mean (SEM) of the S signal
  double meanCerenkov;  // Mean of the total C signal
  double sigmaCerenkov; // Sigma (RMS) of the total C signal
  double semCerenkov;   // Standard Error on the Mean (SEM) of the C signal
  TH1F* h_TotalS;       // Pointer to the total Scintillation histogram
  TH1F* h_TotalC;       // Pointer to the total Cerenkov histogram
};

/**
 * @brief Analyzes a single input file and calculates the total signal (Barrel + Endcap).
 * @param input_file Path to the input PODIO file.
 * @param theta The theta value, used for unique histogram naming.
 * @return An AnalysisResults object with stats and pointers to the created histograms.
 */
AnalysisResults analyzeFile(const std::string& input_file, double theta) {

  std::cout << "  Opening file: " << input_file << std::endl;

  auto reader = podio::makeReader(input_file);

  // Create unique names and titles for histograms
  std::stringstream ss_s_name, ss_c_name, ss_s_title, ss_c_title;
  ss_s_name << "h_TotalS_theta_" << theta;
  ss_c_name << "h_TotalC_theta_" << theta;
  ss_s_title << "Total S (Theta=" << theta << " deg);Total S [p.e.]; Events";
  ss_c_title << "Total C (Theta=" << theta << " deg);Total C [p.e.]; Events";

  // Histograms for the TOTAL signal (Barrel + Endcap)
  // We create them with 'new' so they persist after the function returns
  TH1F* h_TotalS = new TH1F(ss_s_name.str().c_str(), ss_s_title.str().c_str(), 400, 0.0, 12000.0);
  TH1F* h_TotalC = new TH1F(ss_c_name.str().c_str(), ss_c_title.str().c_str(), 300, 0.0, 8000.0);

  // Prevent histograms from being added to any open TFile or TDirectory
  h_TotalS->SetDirectory(nullptr);
  h_TotalC->SetDirectory(nullptr);

  unsigned int nEvents = reader.getEvents();
  std::cout << "  Number of events in file: " << nEvents << std::endl;

  // Loop over events
  for (size_t i = 0; i < nEvents; ++i) {
    if (i > 0 && i % 1000 == 0) {
      std::cout << "    ... processed event " << i << " of " << nEvents << std::endl;
    }

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

    // Calculate the total signal for this event
    float eventTotalS = TotalDRBTS + TotalDRETS;
    float eventTotalC = TotalDRBTC + TotalDRETC;

    // Fill the total histograms
    h_TotalS->Fill(eventTotalS);
    h_TotalC->Fill(eventTotalC);
  } // End event loop

  //reader->closeFile();

  // Calculate mean, sigma (RMS), and SEM (Standard Error of the Mean)
  double nEntriesS = h_TotalS->GetEntries();
  double nEntriesC = h_TotalC->GetEntries();

  AnalysisResults result;
  result.meanScint = h_TotalS->GetMean();
  result.sigmaScint = h_TotalS->GetRMS();
  result.semScint = (nEntriesS > 0) ? (result.sigmaScint / std::sqrt(nEntriesS)) : 0.0; // SEM = RMS / sqrt(N)

  result.meanCerenkov = h_TotalC->GetMean();
  result.sigmaCerenkov = h_TotalC->GetRMS();
  result.semCerenkov = (nEntriesC > 0) ? (result.sigmaCerenkov / std::sqrt(nEntriesC)) : 0.0; // SEM = RMS / sqrt(N)

  // Store pointers to the histograms in the result struct
  result.h_TotalS = h_TotalS;
  result.h_TotalC = h_TotalC;

  return result;
}

/**
 * @brief Main function that orchestrates the analysis over multiple files.
 * This is the function to call in ROOT.
 */
int runHadCalibrationAnalysis() {

  // Map of theta values (in degrees) to their corresponding filenames
  std::map<double, std::string> filesToAnalyze;
  filesToAnalyze[0.5] = "hadcalibration/IDEA_o2_v01_phi0p5_theta0p5.root";
  filesToAnalyze[20.5] = "hadcalibration/IDEA_o2_v01_phi0p5_theta20p5.root";
  filesToAnalyze[40.5] = "hadcalibration/IDEA_o2_v01_phi0p5_theta40p5.root";
  filesToAnalyze[50.6] = "hadcalibration/IDEA_o2_v01_phi0p5_theta50p5.root";
  filesToAnalyze[60.5] = "hadcalibration/IDEA_o2_v01_phi0p5_theta60p5.root";
  filesToAnalyze[70.5] = "hadcalibration/IDEA_o2_v01_phi0p5_theta70p5.root";

  std::vector<double> meanScintValues;
  std::vector<double> meanCerenkovValues;
  std::vector<TH1F*> scintHistograms;    // To store histogram pointers
  std::vector<TH1F*> cerenkovHistograms; // To store histogram pointers

  // --- Initialize TGraphErrors ---
  auto grS = new TGraphErrors();
  grS->SetName("g_MeanS_vs_Theta");
  grS->SetTitle("Mean Scintillation Signal vs. Theta;Theta [degrees];Mean S [p.e.]");
  grS->SetMarkerStyle(20); // Full round marker
  grS->SetMarkerColor(kBlue);
  grS->SetLineColor(kBlue);

  auto grC = new TGraphErrors();
  grC->SetName("g_MeanC_vs_Theta");
  grC->SetTitle("Mean Cerenkov Signal vs. Theta;Theta [degrees];Mean C [p.e.]");
  grC->SetMarkerStyle(20); // Full round marker
  grC->SetMarkerColor(kRed);
  grC->SetLineColor(kRed);

  int pointIndex = 0; // Index for TGraphErrors points

  std::cout << "Starting calibration analysis..." << std::endl;

  for (const auto& pair : filesToAnalyze) {
    double theta = pair.first;
    std::string filename = pair.second;

    // Check if file exists before trying
    if (gSystem->AccessPathName(filename.c_str())) {
      std::cerr << "ERROR: Cannot find file: " << filename << std::endl;
      std::cerr << "Skipping analysis for theta = " << theta << " deg." << std::endl;
      continue; // Skip to the next file
    }

    std::cout << "\nAnalyzing for theta = " << theta << " deg (File: " << filename << ")" << std::endl;

    AnalysisResults results = analyzeFile(filename, theta); // Pass theta

    // Check if the analysis produced valid (non-zero) results
    if (results.meanScint > 0 || results.meanCerenkov > 0) {
      std::cout << "  Results:" << std::endl;
      std::cout << "    -> Mean Scintillation (S): " << results.meanScint << " +/- " << results.semScint
                << " p.e. (SEM)" << std::endl;
      std::cout << "    -> Mean Cerenkov (C):      " << results.meanCerenkov << " +/- " << results.semCerenkov
                << " p.e. (SEM)" << std::endl;

      // Add to vectors for the final average
      meanScintValues.push_back(results.meanScint);
      meanCerenkovValues.push_back(results.meanCerenkov);

      // Store histogram pointers
      if (results.h_TotalS)
        scintHistograms.push_back(results.h_TotalS);
      if (results.h_TotalC)
        cerenkovHistograms.push_back(results.h_TotalC);

      // --- Add to TGraphErrors ---
      // X = theta, Y = mean, X_error = 0, Y_error = SEM
      grS->SetPoint(pointIndex, theta, results.meanScint);
      grS->SetPointError(pointIndex, 0.0, results.semScint); // Null error on X, SEM on Y

      grC->SetPoint(pointIndex, theta, results.meanCerenkov);
      grC->SetPointError(pointIndex, 0.0, results.semCerenkov); // Null error on X, SEM on Y

      pointIndex++; // Increment the point index

    } else {
      std::cerr << "  No valid results from analysis for theta = " << theta << " deg." << std::endl;
      // Clean up histograms if analysis failed
      delete results.h_TotalS;
      delete results.h_TotalC;
    }
  }

  // Calculate the average of the mean values (the "calibration constant")
  double sumMeanS = std::accumulate(meanScintValues.begin(), meanScintValues.end(), 0.0);
  double sumMeanC = std::accumulate(meanCerenkovValues.begin(), meanCerenkovValues.end(), 0.0);

  double averageMeanS = (meanScintValues.size() > 0) ? sumMeanS / meanScintValues.size() : 0.0;
  double averageMeanC = (meanCerenkovValues.size() > 0) ? sumMeanC / meanCerenkovValues.size() : 0.0;

  std::cout << "\n--- Final Result (Average over Theta values) ---" << std::endl;
  if (meanScintValues.size() > 0) {
    std::cout << "Average S signal (based on " << meanScintValues.size() << " points): " << averageMeanS << " p.e."
              << std::endl;
    std::cout << "Average C signal (based on " << meanCerenkovValues.size() << " points): " << averageMeanC << " p.e."
              << std::endl;
  } else {
    std::cout << "No valid files were analyzed. Cannot calculate average." << std::endl;
  }

  // --- Create and save TGraphErrors and Histograms ---
  if (pointIndex > 0) {
    std::cout << "\nSaving TGraphErrors and Histograms to 'CalibrationGraphs.root'..." << std::endl;

    TFile* outFile = new TFile("HadCalibrationGraphs.root", "RECREATE");

    // --- Write Graphs and Canvases ---
    // Draw and save the S graph
    TCanvas* cS = new TCanvas("canvasS", "Scintillation Signal vs. Theta", 800, 600);
    cS->SetGrid();
    grS->Sort();      // Sort points by theta for correct line drawing
    grS->Draw("APL"); // "A" for axes, "P" for markers, "L" for line
    cS->Write();      // Save the TCanvas
    grS->Write();     // Save the TGraphErrors itself

    // Draw and save the C graph
    TCanvas* cC = new TCanvas("canvasC", "Cerenkov Signal vs. Theta", 800, 600);
    cC->SetGrid();
    grC->Sort();      // Sort points by theta
    grC->Draw("APL"); // "A" for axes, "P" for markers, "L" for line
    cC->Write();      // Save the TCanvas
    grC->Write();     // Save the TGraphErrors itself

    // --- Write Histograms ---
    std::cout << "Writing histograms..." << std::endl;
    // Create directories in the ROOT file for organization
    outFile->mkdir("ScintillationHistograms");
    outFile->cd("ScintillationHistograms");
    for (auto h : scintHistograms) {
      if (h)
        h->Write();
    }

    outFile->mkdir("CerenkovHistograms");
    outFile->cd("CerenkovHistograms");
    for (auto h : cerenkovHistograms) {
      if (h)
        h->Write();
    }

    outFile->Close();
    std::cout << "Graphs and histograms saved." << std::endl;

    // --- Clean up memory ---
    // We created these histograms with 'new', so we must delete them
    for (auto h : scintHistograms)
      delete h;
    for (auto h : cerenkovHistograms)
      delete h;
    // TGraphErrors and TCanvas are managed by the TFile,
    // but it's good practice to delete them if they are no longer needed,
    // although ROOT's memory management often handles this.
    // We'll let them be cleaned up by ROOT's gDirectory.
    // Or, more safely:
    delete grS;
    delete grC;
    delete cS;
    delete cC;

  } else {
    std::cout << "\nNo valid points to draw." << std::endl;
    // Clean up empty graphs if they weren't used
    delete grS;
    delete grC;
    // Clean up any histograms that were created but not saved
    for (auto h : scintHistograms)
      delete h;
    for (auto h : cerenkovHistograms)
      delete h;
  }

  return 0;
}
