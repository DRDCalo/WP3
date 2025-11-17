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

// --- RESULTS STRUCTURE ---
// This struct holds the statistics extracted from the Gaussian fit of each histogram.
struct AnalysisResults {
    // Statistics for Calibrated Scintillation Signal (in GeV)
    double meanScint;     // Mean (from gaus fit)
    double sigmaScint;    // Sigma (from gaus fit)
    double semScint;      // Error on Mean (from gaus fit)
    double sigmaScintErr; // Error on Sigma (from gaus fit)
    
    // Statistics for Calibrated Cerenkov Signal (in GeV)
    double meanCerenkov;
    double sigmaCerenkov;
    double semCerenkov; 
    double sigmaCerenkovErr;
    
    // Statistics for Combined Signal (in GeV)
    double meanCombined;
    double sigmaCombined;
    double semCombined;
    double sigmaCombinedErr;
    
    // Pointers to the histograms
    TH1F* h_TotalS;     // Histogram for S [GeV]
    TH1F* h_TotalC;     // Histogram for C [GeV]
    TH1F* h_TotalE;     // Histogram for Combined E [GeV]
};

/**
 * @brief Analyzes a single input file, creates histograms in GeV, and fits them.
 * @param input_file Path to the PODIO file.
 * @param energy The beam energy (GeV), used for histogram naming.
 * @param Ks Scintillation calibration constant (p.e./GeV).
 * @param Kc Cerenkov calibration constant (p.e./GeV).
 * @return An AnalysisResults struct with stats from Gaussian fits and histogram pointers.
 */
AnalysisResults analyzeFile(const std::string& input_file, double energy, double Ks, double Kc) {
    
    std::cout << "  Opening file: " << input_file << std::endl;
    auto reader = podio::makeReader(input_file);

    // --- Create unique names and titles for histograms (all in GeV) ---
    std::stringstream ss_s_name, ss_c_name, ss_e_name, ss_s_title, ss_c_title, ss_e_title;
    ss_s_name << "h_TotalS_E_" << energy;
    ss_c_name << "h_TotalC_E_" << energy;
    ss_e_name << "h_TotalE_E_" << energy;
    
    // Update titles to reflect GeV units
    ss_s_title << "Total S (E=" << energy << " GeV);Total Calibrated S [GeV]; Events";
    ss_c_title << "Total C (E=" << energy << " GeV);Total Calibrated C [GeV]; Events";
    ss_e_title << "Total E Combined (E=" << energy << " GeV);E = (S_{cal} + C_{cal})/2 [GeV]; Events";

    // Create histograms with 'new' and SetDirectory(nullptr) to make them persist.
    // Update ranges to be appropriate for GeV.
    double range_min = energy * 0.5;
    double range_max = energy * 1.5;
    TH1F* h_TotalS = new TH1F(ss_s_name.str().c_str(), ss_s_title.str().c_str(), 200, range_min, range_max); 
    TH1F* h_TotalC = new TH1F(ss_c_name.str().c_str(), ss_c_title.str().c_str(), 200, range_min, range_max);
    TH1F* h_TotalE = new TH1F(ss_e_name.str().c_str(), ss_e_title.str().c_str(), 200, range_min, range_max);

    h_TotalS->SetDirectory(nullptr);
    h_TotalC->SetDirectory(nullptr);
    h_TotalE->SetDirectory(nullptr);

    unsigned int nEvents = reader.getEvents();
    std::cout << "  Number of events in file: " << nEvents << std::endl;

    // --- Event Loop ---
    for (size_t i = 0; i < nEvents; ++i) {
        
        auto event = reader.readNextEvent();

        auto& crystals_Scounts = event.get<edm4hep::SimCalorimeterHitCollection>("SCEPCal_MainScounts");
        auto& crystals_Ccounts = event.get<edm4hep::SimCalorimeterHitCollection>("SCEPCal_MainCcounts");

        float TotalS_pe = 0.;
        float TotalC_pe = 0.;

        for (const auto& crystalS_hit : crystals_Scounts) { TotalS_pe += crystalS_hit.getEnergy(); }
        for (const auto& crystalC_hit : crystals_Ccounts) { TotalC_pe += crystalC_hit.getEnergy(); }
        
        // Calibrated signals per event (in GeV)
        float eventTotalS_cal = TotalS_pe / Ks;
        float eventTotalC_cal = TotalC_pe / Kc;
        
        // Combined Energy (in GeV)
        float eventTotalE = (eventTotalS_cal + eventTotalC_cal) / 2.0;

        // Fill histograms with calibrated GeV values
        h_TotalS->Fill(eventTotalS_cal);
        h_TotalC->Fill(eventTotalC_cal);
        h_TotalE->Fill(eventTotalE);
    } // --- End Event Loop ---
    
    // --- Calculate Statistics from Gaussian Fit ---
    AnalysisResults result;

    // Helper function to fit a histogram and extract parameters
    auto fitHistogram = [&](TH1F* h) -> std::tuple<double, double, double, double> {
        if (h->GetEntries() == 0) {
            std::cerr << "  WARNING: Histogram " << h->GetName() << " is empty. Skipping fit." << std::endl;
            return {0, 0, 0, 0};
        }
        
        // Initial parameter estimates from histogram
        double initialMean = h->GetMean();
        double initialSigma = h->GetRMS();
        double fitMin = initialMean - 2.0 * initialSigma;
        double fitMax = initialMean + 2.0 * initialSigma;

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

    // Fit S histogram
    auto [meanS, sigmaS, semS, sigmaSErr] = fitHistogram(h_TotalS);
    result.meanScint = meanS;
    result.sigmaScint = sigmaS;
    result.semScint = semS;
    result.sigmaScintErr = sigmaSErr;

    // Fit C histogram
    auto [meanC, sigmaC, semC, sigmaCErr] = fitHistogram(h_TotalC);
    result.meanCerenkov = meanC;
    result.sigmaCerenkov = sigmaC;
    result.semCerenkov = semC;
    result.sigmaCerenkovErr = sigmaCErr;

    // Fit E (Combined) histogram
    auto [meanE, sigmaE, semE, sigmaEErr] = fitHistogram(h_TotalE);
    result.meanCombined = meanE;
    result.sigmaCombined = sigmaE;
    result.semCombined = semE;
    result.sigmaCombinedErr = sigmaEErr;

    // Store histogram pointers
    result.h_TotalS = h_TotalS;
    result.h_TotalC = h_TotalC;
    result.h_TotalE = h_TotalE;

    return result;
}

/**
 * @brief Helper function to create a linearity canvas with a ratio plot.
 * @param canvasName Name for the TCanvas.
 * @param canvasTitle Title for the TCanvas.
 * @param graph The TGraphErrors with the linearity points.
 * @param fitFunc The TF1 (pol1) fit to the data.
 * @return A new TGraphErrors object containing the ratio points (caller must delete).
 */
TGraphErrors* createLinearityCanvas(const char* canvasName, const char* canvasTitle, 
                                    TGraphErrors* graph, TF1* fitFunc) 
{
    TCanvas *canvas = new TCanvas(canvasName, canvasTitle, 800, 700);
    
    // Create ratio graph
    TGraphErrors* ratioGraph = new TGraphErrors(graph->GetN());
    ratioGraph->SetName(Form("%s_ratio", graph->GetName()));
    ratioGraph->SetTitle(""); // No title for ratio
    ratioGraph->SetMarkerStyle(graph->GetMarkerStyle());
    ratioGraph->SetMarkerColor(graph->GetMarkerColor());
    ratioGraph->SetLineColor(graph->GetLineColor());

    double x, y, y_fit, y_err, perc_dev, perc_err;
    for (int i = 0; i < graph->GetN(); ++i) {
        graph->GetPoint(i, x, y);
        y_err = graph->GetErrorY(i);
        y_fit = fitFunc->Eval(x);
        
        if (y_fit != 0) {
            perc_dev = (y - y_fit) / y_fit * 100.0;
            perc_err = (y_err / y_fit) * 100.0; // Propagated error
        } else {
            perc_dev = 0;
            perc_err = 0;
        }
        ratioGraph->SetPoint(i, x, perc_dev);
        ratioGraph->SetPointError(i, 0, perc_err);
    }

    // --- Create Pads ---
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0.02); // No bottom margin
    pad1->SetGrid();
    pad1->Draw();
    
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.3);
    pad2->SetTopMargin(0.02); // No top margin
    pad2->SetBottomMargin(0.3); // Make space for X axis label
    pad2->SetGrid();
    pad2->Draw();

    // --- Draw main plot ---
    pad1->cd();
    graph->Draw("AP"); // "AP" = Axes, Points (no line connecting points)
    fitFunc->Draw("SAME");
    graph->GetYaxis()->SetTitleSize(0.05);
    graph->GetYaxis()->SetTitleOffset(0.9);
    graph->GetYaxis()->SetLabelSize(0.04);
    graph->GetXaxis()->SetLabelSize(0); // Hide X label on top plot

    // --- Draw ratio plot ---
    pad2->cd();
    ratioGraph->Draw("AP");
    
    // Style ratio plot axes
    ratioGraph->GetYaxis()->SetTitle("Deviation [%]");
    ratioGraph->GetYaxis()->SetNdivisions(505); // 5 primary, 5 secondary
    ratioGraph->GetYaxis()->SetTitleSize(0.12);
    ratioGraph->GetYaxis()->SetTitleOffset(0.35);
    ratioGraph->GetYaxis()->SetLabelSize(0.1);
    ratioGraph->GetYaxis()->CenterTitle();
    
    ratioGraph->GetXaxis()->SetTitle(graph->GetXaxis()->GetTitle()); // Get title from main graph
    ratioGraph->GetXaxis()->SetTitleSize(0.12);
    ratioGraph->GetXaxis()->SetTitleOffset(1.0);
    ratioGraph->GetXaxis()->SetLabelSize(0.1);

    // Set Y range for ratio, e.g., +/- 5%
    ratioGraph->SetMinimum(-1.0);
    ratioGraph->SetMaximum(1.0);

    // Draw line at 0
    TLine *line = new TLine(pad2->GetUxmin(), 0, pad2->GetUxmax(), 0);
    line->SetLineColor(kRed);
    line->SetLineStyle(2); // Dashed
    line->Draw();
    
    canvas->Write();
    ratioGraph->Write();
    
    // Return ratio graph so it can be deleted later
    return ratioGraph;
}


/**
 * @brief Main function to run the analysis and produce linearity and resolution graphs.
 */
int emResolutionFit() {
    
    // --- CALIBRATION CONSTANTS (p.e./GeV) ---
    const double CALIBRATION_CONSTANT_S = 1960.92; // S p.e./GeV
    const double CALIBRATION_CONSTANT_C = 97.531; // C p.e./GeV
    
    // --- FILE MAP (Energy in GeV) ---
    std::map<double, std::string> filesToAnalyze;
    filesToAnalyze[10.0] = "emlinearity/IDEA_o2_v01_phi0p5_theta0p5_10gev.root";
    filesToAnalyze[20.0] = "emlinearity/IDEA_o2_v01_phi0p5_theta0p5_20gev.root";
    filesToAnalyze[40.0] = "emlinearity/IDEA_o2_v01_phi0p5_theta0p5_40gev.root";
    // ... add more energies and files here

    std::vector<TH1F*> allHistograms;
    std::vector<TObject*> objectsToClean; // Store all created ROOT objects

    // --- 1. Linearity Graphs (Mean Signal [GeV] vs. True Energy [GeV]) ---
    TGraphErrors* grS_Linearity = new TGraphErrors();
    grS_Linearity->SetName("g_Linearity_S");
    grS_Linearity->SetTitle("Scintillation Signal Linearity;True Energy [GeV];Mean S [GeV]");
    grS_Linearity->SetMarkerStyle(20); grS_Linearity->SetMarkerColor(kBlue); grS_Linearity->SetLineColor(kBlue);
    objectsToClean.push_back(grS_Linearity);

    TGraphErrors* grC_Linearity = new TGraphErrors();
    grC_Linearity->SetName("g_Linearity_C");
    grC_Linearity->SetTitle("Cerenkov Signal Linearity;True Energy [GeV];Mean C [GeV]");
    grC_Linearity->SetMarkerStyle(20); grC_Linearity->SetMarkerColor(kRed); grC_Linearity->SetLineColor(kRed);
    objectsToClean.push_back(grC_Linearity);

    TGraphErrors* grE_Linearity = new TGraphErrors();
    grE_Linearity->SetName("g_Linearity_Combined");
    grE_Linearity->SetTitle("Combined Signal Linearity;True Energy [GeV];Mean E=(S_{cal}+C_{cal})/2 [GeV]");
    grE_Linearity->SetMarkerStyle(20); grE_Linearity->SetMarkerColor(kGreen+2); grE_Linearity->SetLineColor(kGreen+2);
    objectsToClean.push_back(grE_Linearity);

    // --- 2. Resolution Graphs (Sigma/Mean vs. Energy) ---
    TGraphErrors* grS_Resolution = new TGraphErrors();
    grS_Resolution->SetName("g_Resolution_S");
    grS_Resolution->SetTitle("Scintillation Energy Resolution;Energy [GeV];Resolution (#sigma/#mu)");
    grS_Resolution->SetMarkerStyle(20); grS_Resolution->SetMarkerColor(kBlue); grS_Resolution->SetLineColor(kBlue);
    objectsToClean.push_back(grS_Resolution);

    TGraphErrors* grC_Resolution = new TGraphErrors();
    grC_Resolution->SetName("g_Resolution_C");
    grC_Resolution->SetTitle("Cerenkov Energy Resolution;Energy [GeV];Resolution (#sigma/#mu)");
    grC_Resolution->SetMarkerStyle(20); grC_Resolution->SetMarkerColor(kRed); grC_Resolution->SetLineColor(kRed);
    objectsToClean.push_back(grC_Resolution);

    TGraphErrors* grE_Resolution = new TGraphErrors();
    grE_Resolution->SetName("g_Resolution_Combined");
    grE_Resolution->SetTitle("Combined Energy Resolution;Energy [GeV];Resolution (#sigma_{E}/#mu_{E})");
    grE_Resolution->SetMarkerStyle(20); grE_Resolution->SetMarkerColor(kGreen+2); grE_Resolution->SetLineColor(kGreen+2);
    objectsToClean.push_back(grE_Resolution);


    int pointIndex = 0; 
    double minEnergy = 1e9, maxEnergy = -1;

    std::cout << "Starting energy scan analysis with calibration..." << std::endl;
    std::cout << "Calibration S: " << CALIBRATION_CONSTANT_S << " p.e./GeV" << std::endl;
    std::cout << "Calibration C: " << CALIBRATION_CONSTANT_C << " p.e./GeV" << std::endl;

    for (const auto& pair : filesToAnalyze) {
        double energy = pair.first;
        std::string filename = pair.second;

        if (energy < minEnergy) minEnergy = energy;
        if (energy > maxEnergy) maxEnergy = energy;

        if (gSystem->AccessPathName(filename.c_str())) {
            std::cerr << "ERROR: Cannot find file: " << filename << std::endl;
            continue; 
        }

        std::cout << "\nAnalyzing for energy = " << energy << " GeV (File: " << filename << ")" << std::endl;
        
        AnalysisResults results = analyzeFile(filename, energy, CALIBRATION_CONSTANT_S, CALIBRATION_CONSTANT_C);

        if (results.meanScint > 0 && results.meanCerenkov > 0) {
            
            double resS = results.sigmaScint / results.meanScint;
            double resC = results.sigmaCerenkov / results.meanCerenkov;
            double resE = results.sigmaCombined / results.meanCombined;

            // --- Calculate resolution error ---
            // R = sigma / mu
            // (delta_R/R)^2 = (delta_sigma/sigma)^2 + (delta_mu/mu)^2
            auto calc_res_err = [](double res, double sigma, double sigma_err, double mu, double mu_err) {
                if (mu == 0.0 || sigma == 0.0) return 0.0;
                double rel_err_sigma_sq = std::pow(sigma_err / sigma, 2);
                double rel_err_mu_sq = std::pow(mu_err / mu, 2);
                return res * std::sqrt(rel_err_sigma_sq + rel_err_mu_sq);
            };
            
            double err_resS = calc_res_err(resS, results.sigmaScint, results.sigmaScintErr, results.meanScint, results.semScint);
            double err_resC = calc_res_err(resC, results.sigmaCerenkov, results.sigmaCerenkovErr, results.meanCerenkov, results.semCerenkov);
            double err_resE = calc_res_err(resE, results.sigmaCombined, results.sigmaCombinedErr, results.meanCombined, results.semCombined);
            
            // --- 1. Linearity Points ---
            grS_Linearity->SetPoint(pointIndex, energy, results.meanScint);
            grS_Linearity->SetPointError(pointIndex, 0.0, results.semScint);

            grC_Linearity->SetPoint(pointIndex, energy, results.meanCerenkov);
            grC_Linearity->SetPointError(pointIndex, 0.0, results.semCerenkov);
            
            grE_Linearity->SetPoint(pointIndex, energy, results.meanCombined);
            grE_Linearity->SetPointError(pointIndex, 0.0, results.semCombined);

            // --- 2. Resolution Points ---
            grS_Resolution->SetPoint(pointIndex, energy, resS);
            grS_Resolution->SetPointError(pointIndex, 0.0, err_resS);

            grC_Resolution->SetPoint(pointIndex, energy, resC);
            grC_Resolution->SetPointError(pointIndex, 0.0, err_resC);
            
            grE_Resolution->SetPoint(pointIndex, energy, resE);
            grE_Resolution->SetPointError(pointIndex, 0.0, err_resE);
            
            pointIndex++; 

            // Store histogram pointers for saving
            if (results.h_TotalS) allHistograms.push_back(results.h_TotalS);
            if (results.h_TotalC) allHistograms.push_back(results.h_TotalC);
            if (results.h_TotalE) allHistograms.push_back(results.h_TotalE);

        } else {
            std::cerr << "  No valid results from analysis for energy = " << energy << " GeV." << std::endl;
            delete results.h_TotalS;
            delete results.h_TotalC;
            delete results.h_TotalE;
        }
    }

    // --- Create and save TGraphErrors, Fits, and Histograms ---
    if (pointIndex > 0) {
        std::cout << "\nSaving Calibrated Graphs, Fits, and Histograms to 'EnergyScanCalibratedResults.root'..." << std::endl;
        
        gStyle->SetOptFit(0); // Turn off global fit stats box
        TFile *outFile = new TFile("EnergyScanCalibratedResults.root", "RECREATE");

        // Set fit range
        minEnergy = filesToAnalyze.begin()->first * 0.9;
        maxEnergy = filesToAnalyze.rbegin()->first * 1.1;

        // --- Linearity ---
        outFile->mkdir("Linearity");
        outFile->cd("Linearity");

        // FIT (pol1: p0 + p1*x)
        TF1* fit_lin_S = new TF1("fit_lin_S", "pol1", minEnergy, maxEnergy);
        fit_lin_S->SetLineColor(kBlue);
        grS_Linearity->Fit(fit_lin_S, "RQ"); // Fit quietly
        objectsToClean.push_back(fit_lin_S);
        
        TF1* fit_lin_C = new TF1("fit_lin_C", "pol1", minEnergy, maxEnergy);
        fit_lin_C->SetLineColor(kRed);
        grC_Linearity->Fit(fit_lin_C, "RQ"); // Fit quietly
        objectsToClean.push_back(fit_lin_C);
        
        TF1* fit_lin_E = new TF1("fit_lin_E", "pol1", minEnergy, maxEnergy);
        fit_lin_E->SetLineColor(kGreen+2);
        grE_Linearity->Fit(fit_lin_E, "RQ"); // Fit quietly
        objectsToClean.push_back(fit_lin_E);
        
        // Create ratio plots
        TGraphErrors* grS_Ratio = createLinearityCanvas("c_Linearity_S", "Scintillation Linearity", grS_Linearity, fit_lin_S);
        TGraphErrors* grC_Ratio = createLinearityCanvas("c_Linearity_C", "Cerenkov Linearity", grC_Linearity, fit_lin_C);
        TGraphErrors* grE_Ratio = createLinearityCanvas("c_Linearity_Combined", "Combined Energy Linearity", grE_Linearity, fit_lin_E);
        objectsToClean.push_back(grS_Ratio);
        objectsToClean.push_back(grC_Ratio);
        objectsToClean.push_back(grE_Ratio);
        
        grS_Linearity->Write();
        grC_Linearity->Write();
        grE_Linearity->Write();

        // --- Resolution ---
        outFile->mkdir("Resolution");
        outFile->cd("Resolution");
        
        // FIT (Resolution: [0]/sqrt(E) + [1])
        TF1* fit_res_S = new TF1("fit_res_S", "[0]/sqrt(x) + [1]", minEnergy, maxEnergy);
        fit_res_S->SetParNames("Stochastic", "Constant");
        fit_res_S->SetParameters(0.1, 0.01); 
        fit_res_S->SetLineColor(kBlue);
        grS_Resolution->Fit(fit_res_S, "RQ");
        objectsToClean.push_back(fit_res_S);

        TF1* fit_res_C = new TF1("fit_res_C", "[0]/sqrt(x) + [1]", minEnergy, maxEnergy);
        fit_res_C->SetParNames("Stochastic", "Constant");
        fit_res_C->SetParameters(0.1, 0.01); 
        fit_res_C->SetLineColor(kRed);
        grC_Resolution->Fit(fit_res_C, "RQ");
        objectsToClean.push_back(fit_res_C);

        TF1* fit_res_E = new TF1("fit_res_E", "[0]/sqrt(x) + [1]", minEnergy, maxEnergy);
        fit_res_E->SetParNames("Stochastic", "Constant");
        fit_res_E->SetParameters(0.1, 0.01); 
        fit_res_E->SetLineColor(kGreen+2);
        grE_Resolution->Fit(fit_res_E, "RQ");
        objectsToClean.push_back(fit_res_E);

        // S Resolution Canvas
        TCanvas *cS_res = new TCanvas("c_Resolution_S", "Scintillation Resolution", 800, 600);
        cS_res->SetLeftMargin(0.15); // Fix Y label cutoff
        cS_res->SetGrid();
        grS_Resolution->Draw("AP"); // "AP" = Axes, Points
        fit_res_S->Draw("SAME");
        cS_res->Write();
        grS_Resolution->Write();
        objectsToClean.push_back(cS_res);

        // C Resolution Canvas
        TCanvas *cC_res = new TCanvas("c_Resolution_C", "Cerenkov Resolution", 800, 600);
        cC_res->SetLeftMargin(0.15); // Fix Y label cutoff
        cC_res->SetGrid();
        grC_Resolution->Draw("AP");
        fit_res_C->Draw("SAME");
        cC_res->Write();
        grC_Resolution->Write();
        objectsToClean.push_back(cC_res);
        
        // E Combined Resolution Canvas
        TCanvas *cE_res = new TCanvas("c_Resolution_Combined", "Combined Energy Resolution", 800, 600);
        cE_res->SetLeftMargin(0.15); // Fix Y label cutoff
        cE_res->SetGrid();
        grE_Resolution->Draw("AP");
        fit_res_E->Draw("SAME");
        cE_res->Write();
        grE_Resolution->Write();
        objectsToClean.push_back(cE_res);
        
	// --- NEW: S and C Combined Resolution Plot ---
        TCanvas *cSC_res = new TCanvas("c_Resolution_SC_Combined", "S and C Resolution", 800, 600);
        cSC_res->SetLeftMargin(0.15); // Fix Y label cutoff
        cSC_res->SetGrid();
        
        // Use a TMultiGraph to auto-fit the Y-axis range for both graphs
        TMultiGraph *mg_res = new TMultiGraph();
        mg_res->Add(grS_Resolution, "P"); // "P" = draw points
        mg_res->Add(grC_Resolution, "P");
        objectsToClean.push_back(mg_res); // Add to cleanup list

        // Draw the multigraph first. "A" draws axes.
        mg_res->Draw("A"); 
        
        // Set titles on the TMultiGraph's axes
        mg_res->SetTitle("S and C Energy Resolution;Energy [GeV];Resolution (#sigma/#mu)");
        mg_res->GetYaxis()->SetTitleOffset(1.25); // Adjust offset after margin change

        // Now draw the fits on top
        fit_res_S->Draw("SAME");
        fit_res_C->Draw("SAME");

        // Redraw points on top of grid/fits (optional but looks better)
        grS_Resolution->Draw("P SAME");
        grC_Resolution->Draw("P SAME");

        // Legend
        TLegend *leg = new TLegend(0.55, 0.75, 0.88, 0.88);
        leg->AddEntry(grS_Resolution, "Scintillation (S)", "p");
        leg->AddEntry(grC_Resolution, "Cerenkov (C)", "p");
        leg->SetBorderSize(1);
        leg->Draw();
        objectsToClean.push_back(leg);

        // PaveText for S Fit Results
        TPaveText *ptS = new TPaveText(0.18, 0.75, 0.53, 0.88, "NDC");
        ptS->SetBorderSize(1);
        ptS->SetFillColor(0);
        ptS->SetTextColor(kBlue);
        ptS->SetTextAlign(12); // Left-aligned
        ptS->AddText(Form("S Fit: a = %.3f #pm %.4f", fit_res_S->GetParameter(0), fit_res_S->GetParError(0)));
        ptS->AddText(Form("S Fit: c = %.3f #pm %.4f", fit_res_S->GetParameter(1), fit_res_S->GetParError(1)));
        ptS->Draw();
        objectsToClean.push_back(ptS);

        // PaveText for C Fit Results
        TPaveText *ptC = new TPaveText(0.18, 0.60, 0.53, 0.73, "NDC");
        ptC->SetBorderSize(1);
        ptC->SetFillColor(0);
        ptC->SetTextColor(kRed);
        ptC->SetTextAlign(12); // Left-aligned
        ptC->AddText(Form("C Fit: a = %.3f #pm %.4f", fit_res_C->GetParameter(0), fit_res_C->GetParError(0)));
        ptC->AddText(Form("C Fit: c = %.3f #pm %.4f", fit_res_C->GetParameter(1), fit_res_C->GetParError(1)));
        ptC->Draw();
        objectsToClean.push_back(ptC);

        cSC_res->Write();
        objectsToClean.push_back(cSC_res);
        

        // --- Write Histograms ---
        // All histograms are now in GeV
        outFile->mkdir("Histograms_GeV");
        outFile->cd("Histograms_GeV");
        for (auto h : allHistograms) { if (h) h->Write(); }

        outFile->Close();
        std::cout << "Analysis results saved to 'EnergyScanCalibratedResults.root'." << std::endl;
        
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
