#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <mpi.h>
#include <fstream>
#include "goldenmodesB.h"
#include "CKM.h"
#include <vector>
#include <random>
#include <chrono>
#include <thread>
#include <string>
#include <cstring>

using namespace std;

int main(int argc, char* argv[]) {
    // Parse command line arguments
    string outputFileName = "";
    bool flagBJPSIV = false;
    bool flagBJPSIP = false;
    bool flagBDDb = false;
    double dsu3_limit = 0.2;
    double ewp_limit = 0.0;
    int nIterationsPreRun = 100000;
    int nIterationsRun = 10000;
    int nIterationsPreRunFactorized = 1000;
    int nIterationsUpdateMax = 100;

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--outfile") == 0 && i + 1 < argc) {
            outputFileName = argv[++i];
        } else if (strcmp(argv[i], "--flagBJPSIV") == 0) {
            flagBJPSIV = true;
        } else if (strcmp(argv[i], "--flagBJPSIP") == 0) {
            flagBJPSIP = true;
        } else if (strcmp(argv[i], "--flagBDDb") == 0) {
            flagBDDb = true;
        } else if (strcmp(argv[i], "--dsu3_limit") == 0 && i + 1 < argc) {
            dsu3_limit = atof(argv[++i]);
        } else if (strcmp(argv[i], "--ewp_limit") == 0 && i + 1 < argc) {
            ewp_limit = atof(argv[++i]);
        } else if (strcmp(argv[i], "--nPreRun") == 0 && i + 1 < argc) {
            nIterationsPreRun = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--nRun") == 0 && i + 1 < argc) {
            nIterationsRun = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--nPreRunFactorized") == 0 && i + 1 < argc) {
            nIterationsPreRunFactorized = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--nUpdateMax") == 0 && i + 1 < argc) {
            nIterationsUpdateMax = atoi(argv[++i]);
        }
    }

    cout << "Output file name: " << outputFileName << endl;
    cout << "FlagBJPSIV: " << (flagBJPSIV ? "true" : "false") << endl;
    cout << "FlagBJPSIP: " << (flagBJPSIP ? "true" : "false") << endl;
    cout << "FlagBDDb: " << (flagBDDb ? "true" : "false") << endl;
    cout << "dsu3_limit: " << dsu3_limit << endl;
    cout << "ewp_limit: " << ewp_limit << endl;
    cout << "nIterationsPreRun: " << nIterationsPreRun << endl;
    cout << "nIterationsRun: " << nIterationsRun << endl;
    cout << "nIterationsPreRunFactorized: " << nIterationsPreRunFactorized << endl;
    cout << "nIterationsUpdateMax: " << nIterationsUpdateMax << endl;
    // Initialize MPI
    MPI_Init(NULL, NULL);
    // Initialize BAT logging
    BCLog::OpenLog(outputFileName+"log.txt", BCLog::detail, BCLog::detail);
    BCLog::SetLogLevelScreen(BCLog::summary);

    // Create an instance of the BqDqDqbar model
    goldenmodesB model(dsu3_limit, ewp_limit, flagBJPSIP, flagBJPSIV, flagBDDb);
    cout << "Model constructed correctly" << endl;

    // Set the number of chains, pre-run, and run iterations
    model.SetNChains(12);
    model.SetNIterationsPreRunMax(nIterationsPreRun); // Pre-run iterations
    model.SetNIterationsRun(nIterationsRun);
    model.SetProposeMultivariate(false);
    //model.SetNIterationsPreRunFactorized(nIterationsPreRunFactorized);
    model.SetNIterationsPreRunCheck(nIterationsUpdateMax);
    
    // Start time measurement
    auto start = chrono::high_resolution_clock::now();

    // Show parameter ranges
    for (unsigned i = 0; i < model.GetNParameters(); ++i) {
        cout << "Parameter " << i << ": " << model.GetParameter(i).GetName()
             << " with range [" << model.GetParameter(i).GetLowerLimit()
             << ", " << model.GetParameter(i).GetUpperLimit() << "]" << endl;
    }


    model.MarginalizeAll();  // MCMC

    model.PrintSummary();
    model.PrintAllMarginalized(outputFileName+"marginalized_parameters.pdf");
    model.PrintCorrelationMatrix(outputFileName+"correlation_matrix.pdf");
    model.WriteMarginalizedDistributions(outputFileName+"marginalized_pars.root", "RECREATE");
    model.SaveHistograms(outputFileName+"output_histograms.root");
    // model.PrintObservablePulls("observable_pulls_GM.txt");
    // model.PrintHistogram();

    BCLog::OutSummary("MCMC analysis completed.");
    BCLog::CloseLog();

    // End time measurement
    auto endTime = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsedTime = endTime - start;
    cout << "\nExecution time: " << elapsedTime.count() << " seconds" << endl;

    MPI_Finalize();
    return 0;
}
