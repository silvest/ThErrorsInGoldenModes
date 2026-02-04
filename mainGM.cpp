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
    bool Univariate = false;

    double dsu3_limit = 0.2;
    double ewp_limit = 0.0;
    int nIterationsPreRun = 100000;
    int nIterationsRun = 10000;
    int nIterationsPreRunFactorized = 10000;
    int nIterationsUpdateMax = 100;
    int nChains = 10;

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
        else if (strcmp(argv[i], "--nChains") == 0 && i + 1 < argc) {
             nChains = atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "--Univariate") == 0) {
            Univariate = true;
        }
        else if (strcmp(argv[i], "--help") == 0) {
            cout << "Usage: " << argv[0] << " [options]\n"
                 << "Options:\n"
                 << "  --outfile <filename>            Specify output file name prefix\n"
                 << "  --flagBJPSIV                    Enable BJPSIV flag\n"
                 << "  --flagBJPSIP                    Enable BJPSIP flag\n"
                 << "  --flagBDDb                     Enable BDDb flag\n"
                 << "  --dsu3_limit <value>            Set dsu3 limit (default: 0.2)\n"
                 << "  --ewp_limit <value>             Set ewp limit (default: 0.0)\n"
                 << "  --nPreRun <value>               Set number of pre-run iterations (default: 100000)\n"
                 << "  --nRun <value>                  Set number of run iterations (default: 10000)\n"
                 << "  --nPreRunFactorized <value>     Set number of factorized pre-run iterations (default: 1000)\n"
                  "  --nUpdateMax value            Set max number of update iterations (default: 100)\n"
                  "  --nChains value               Set number of chains (default: 10)\n"
                  "  --Univariate                   Enable univariate mode\n"
                  "  --help                          Show this help message\n";
            return 0;
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
    model.SetNChains(nChains);
    model.SetNIterationsPreRunMax(nIterationsPreRun); // Pre-run iterations
    model.SetNIterationsRun(nIterationsRun);
    if (Univariate) {
        model.SetProposeMultivariate(false);
    } else {
        model.SetProposeMultivariate(true);
    }
    model.SetNIterationsPreRunFactorized(nIterationsPreRunFactorized);
    model.SetNIterationsPreRunCheck(nIterationsUpdateMax);
    model.SetInitialPositionScheme(BCEngineMCMC::kInitCenter);
//    model.WriteMarkovChain(outputFileName+"markov_chain.root", "RECREATE", true, true);

    
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
