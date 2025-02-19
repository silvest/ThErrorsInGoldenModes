#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <fstream>
#include "goldenmodes.h"
#include "CKM.h"
#include <vector>
#include <random>
#include <chrono>
#include <thread>

using namespace std;

// Function to show the progress bar
void showProgressBar(int current, int total, double elapsedTime) {
    static int lastPrintedPercentage = -1;  // Store the last printed percentage
    int percentage = (double)current / total * 100;

    // Only update every 10% progress
    if (percentage / 20 > lastPrintedPercentage / 20) {
        lastPrintedPercentage = percentage;  // Update last printed percentage
        double remainingTime = (elapsedTime / current) * (total - current);

        // Print progress bar and estimated time remaining
        cout << "\rProgress: [" << string(percentage / 2, '=') << string(50 - (percentage / 2), ' ')
             << "] " << percentage << "%  ";
        cout << "Estimated time remaining: " << remainingTime << " seconds" << flush;
    }
}


int main() {
    // Initialize BAT logging
    BCLog::OpenLog("log_GM.txt", BCLog::detail, BCLog::detail);
    BCLog::SetLogLevelScreen(BCLog::summary);

    // Create an instance of the BqDqDqbar model
    goldenmodes model;
    cout << "Model constructed correctly" << endl;

    // Set the number of chains, pre-run, and run iterations
    model.SetNChains(4);
    model.SetNIterationsPreRunMax(2000000); // Pre-run iterations
    model.SetNIterationsRun(1000000);
    model.SetProposeMultivariate(true);
    unsigned int nParameters = model.GetNParameters();
    cout << "Expected number of parameters: " << nParameters << endl;

    model.SetNIterationsPreRunCheck(1000);
    model.SetInitialPositionScheme(BCEngineMCMC::kInitCenter);
    cout << "Set number of chains, pre-runs, and run iterations..." << endl;

    // Start time measurement
    auto start = chrono::high_resolution_clock::now();

    // Show parameter ranges
    for (unsigned i = 0; i < model.GetNParameters(); ++i) {
        cout << "Parameter " << i << ": " << model.GetParameter(i).GetName()
             << " with range [" << model.GetParameter(i).GetLowerLimit()
             << ", " << model.GetParameter(i).GetUpperLimit() << "]" << endl;
    }

    // ==============================
    // PRE-RUN Progress Tracking
    // ==============================
    cout << "Starting Pre-run..." << endl;
    unsigned int totalPreRun = model.GetNIterationsPreRunMax();
    unsigned int updateFrequencyPreRun = totalPreRun / 100; // Update every 1%

    for (unsigned int i = 0; i < totalPreRun; i += updateFrequencyPreRun) {
        model.FindMode(model.GetBestFitParameters()); // Pre-run step

        auto now = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsedTime = now - start;
        showProgressBar(i, totalPreRun, elapsedTime.count());
    }

    showProgressBar(totalPreRun, totalPreRun, 0); // Ensure completion is printed
    cout << endl;
    cout << "Pre-run completed!" << endl;

    // ==============================
    // MCMC RUN Progress Tracking
    // ==============================
    cout << "Starting MCMC Run..." << endl;
    unsigned int totalRun = model.GetNIterationsRun();
    unsigned int updateFrequencyRun = totalRun / 100; // Update every 1%

    for (unsigned int i = 0; i < totalRun; i += updateFrequencyRun) {
        model.MarginalizeAll();  // MCMC step

        auto now = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsedTime = now - start;
        showProgressBar(i, totalRun, elapsedTime.count());
    }

    showProgressBar(totalRun, totalRun, 0); // Ensure completion is printed
    cout << endl;
    cout << "MCMC Run completed!" << endl;

    // ==============================
    // Save Results & Log Output
    // ==============================
    std::vector<double> bestFitParameters = model.GetBestFitParameters();

    std::ofstream resultFile("results_GM.txt");
    resultFile << "Best-fit parameters:\n";
    for (unsigned i = 0; i < model.GetNParameters(); ++i) {
        resultFile << model.GetParameter(i).GetName() << ": "
                   << model.GetBestFitParameters()[i] << std::endl;
    }
    resultFile.close();

    model.PrintSummary();
    model.PrintAllMarginalized("marginalized_parameters_GM.pdf");
    model.WriteMarginalizedDistributions("marginalized_pars_GM.root", "RECREATE");
    model.SaveHistograms("output_histograms_GM.root");
    model.PrintObservablePulls("observable_pulls_GM.txt");

    BCLog::OutSummary("MCMC analysis completed.");
    BCLog::CloseLog();

    // End time measurement
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsedTime = end - start;
    cout << "\nExecution time: " << elapsedTime.count() << " seconds" << endl;

    return 0;
}
