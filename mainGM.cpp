#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <mpi.h>
#include <fstream>
#include "goldenmodesB.h"
#include "goldenmodesB_indSU3.h"
#include "CKM.h"
#include <vector>
#include <random>
#include <chrono>
#include <thread>
#include <string>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <limits>

using namespace std;

// Template helper: run MCMC and save results for any model sharing the BCModel+SaveHistograms interface
template <typename Model>
void runModel(Model &model, const string &outputFileName,
              int nChains, int nIterationsPreRun, int nIterationsRun,
              int nIterationsPreRunFactorized, int nIterationsUpdateMax,
              bool Univariate, const vector<double>& temperatures, unsigned swapInterval,
              const string& proposalInput, const string& proposalOutput, unsigned proposalWarmup,
              bool ptDiagnostics, unsigned diagnosticInterval, unsigned traceInterval)
{
    model.SetNChains(nChains);
    if (!temperatures.empty()) {
        model.SetParallelTempering(temperatures, nChains, swapInterval);
        if (ptDiagnostics) model.SetPTDiagnosticsOutput(outputFileName.empty() ? "./" : outputFileName, diagnosticInterval, traceInterval);
        cout << "Parallel Tempering: " << model.GetNColdChains() << " independent ladders, "
             << temperatures.size() << " temperatures, " << model.GetNChains() << " total replicas" << endl;
        model.WriteMarkovChain(outputFileName + "cold_chains.root", "RECREATE", true, false);
        ofstream config(outputFileName + "pt_temperatures.csv");
        if (!config) throw runtime_error("Cannot write PT temperature configuration");
        config.precision(17);
        config << "slot,ladder,temperature\n";
        for (unsigned c = 0; c < model.GetNChains(); ++c)
            config << c << ',' << c % model.GetNColdChains() << ',' << model.GetChainTemperature(c) << '\n';
    }
    model.SetNIterationsPreRunMax(nIterationsPreRun);
    model.SetNIterationsRun(nIterationsRun);
    model.SetProposeMultivariate(!Univariate);
    model.SetNIterationsPreRunFactorized(nIterationsPreRunFactorized);
    model.SetNIterationsPreRunCheck(nIterationsUpdateMax);
    model.SetInitialPositionScheme(BCEngineMCMC::kInitRandomPrior);
    model.SetProposalFunctionOutput(proposalOutput.empty() ? outputFileName + "proposal.txt" : proposalOutput);
    if (!proposalInput.empty()) model.LoadProposalFunction(proposalInput, proposalWarmup);

    for (unsigned i = 0; i < model.GetNParameters(); ++i)
        cout << "Parameter " << i << ": " << model.GetParameter(i).GetName()
             << " [" << model.GetParameter(i).GetLowerLimit()
             << ", " << model.GetParameter(i).GetUpperLimit() << "]" << endl;

    if (!Univariate && nIterationsPreRunFactorized >= nIterationsPreRun && proposalInput.empty())
        BCLog::OutWarning("Factorized warm-up covers the entire prerun: the joint multivariate proposal will have no prerun sampling phase.");
    model.MarginalizeAll();
    if (!temperatures.empty()) {
        model.WriteParallelTemperingDiagnostics(outputFileName + "pt_swaps.csv");
        model.WriteParallelTemperingTransport(outputFileName);
        model.CloseOutputFile();
    }
    model.FindMode();
    model.PrintSummary();
    model.PrintAllMarginalized(outputFileName + "marginalized_parameters.pdf");
    model.PrintCorrelationMatrix(outputFileName + "correlation_matrix.pdf");
    model.WriteMarginalizedDistributions(outputFileName + "marginalized_pars.root", "RECREATE");
    model.SaveHistograms(outputFileName + "output_histograms.root");
}

// Worker loop: non-rank-0 processes call this and compute LogLikelihood on demand.
template <typename Model>
void workerLoop(Model &model)
{
    const int buffsize = static_cast<int>(model.GetNParameters()) + 1;
    vector<double> recvbuff(buffsize);
    vector<double> pars(model.GetNParameters());
    double ll;
    while (true) {
        MPI_Scatter(nullptr, buffsize, MPI_DOUBLE,
                    recvbuff.data(), buffsize, MPI_DOUBLE,
                    0, MPI_COMM_WORLD);
        if (recvbuff[0] == 1.) {          // evaluate log-likelihood
            pars.assign(recvbuff.begin() + 1, recvbuff.end());
            ll = model.LogEval(pars);
            MPI_Gather(&ll, 1, MPI_DOUBLE, nullptr, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        } else if (recvbuff[0] == 0.) {   // dummy slot — return -inf
            ll = log(0.);
            MPI_Gather(&ll, 1, MPI_DOUBLE, nullptr, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        } else if (recvbuff[0] == -1.) {  // termination
            break;
        }
    }
}

// Send termination signal to all worker processes after MarginalizeAll.
void terminateWorkers(int procnum, int buffsize)
{
    if (procnum <= 1) return;
    vector<double> sendbuff(static_cast<size_t>(procnum) * buffsize, 0.);
    for (int il = 1; il < procnum; ++il)
        sendbuff[il * buffsize] = -1.;   // task code -1 = exit
    vector<double> dummy(buffsize);
    MPI_Scatter(sendbuff.data(), buffsize, MPI_DOUBLE,
                dummy.data(), buffsize, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
}

int main(int argc, char* argv[]) {
    // Parse command line arguments
    string outputFileName = "";
    bool flagBJPSIV = false;
    bool flagBJPSIP = false;
    bool flagBDDb = false;
    bool Univariate = false;
    bool flagIndSU3 = false;
    bool ptFullPriors = false;
    bool flagSU3ReIm = false;
    bool flagGaussianCKM = false;
    bool flagPositiveG2tP = false;
    bool flagPositiveEA1P = false;
    bool flagPositiveImEA2tP = false;
    bool flagSmallE2tP = false;
    bool flagMinimal = false;
    bool flagPositiveE2tV = false;
    bool flagPositiveG2tV = false;
    bool flagSmallEA1P = false;
    bool flagPositiveG1tD = false;
    bool flagPositiveA1D = false;
    string su3_weight = "abs"; // default SU(3) weight mode for indSU3

    double dsu3_limit = 0.2;
    double ewp_limit = 0.0;
    double su3_sigma = -1.0;  // <=0 means free BAT parameter
    int nIterationsPreRun = 100000;
    int nIterationsRun = 10000;
    int nIterationsPreRunFactorized = 10000;
    int nIterationsUpdateMax = 100;
    int nChains = 10;
    vector<double> temperatures;
    unsigned swapInterval = 10;
    bool ptDiagnostics = true, diagnosticOptions = false;
    unsigned diagnosticInterval = 1000, traceInterval = 100;
    string proposalInput, proposalOutput;
    unsigned proposalWarmup = 1000;

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--no-pt-diagnostics") == 0) {
            ptDiagnostics = false; diagnosticOptions = true;
        } else if (strcmp(argv[i], "--pt-diagnostics-interval") == 0 || strcmp(argv[i], "--pt-trace-interval") == 0) {
            const bool trace = strcmp(argv[i], "--pt-trace-interval") == 0;
            if (i+1 >= argc) { cerr << argv[i] << " requires an integer" << endl; return 1; }
            const string value = argv[++i];
            try {
                size_t used = 0;
                const long long n = stoll(value, &used);
                if (used != value.size() || n < (trace ? 0 : 1) || n > numeric_limits<unsigned>::max())
                    throw invalid_argument("out of range");
                (trace ? traceInterval : diagnosticInterval) = static_cast<unsigned>(n);
                diagnosticOptions = true;
            } catch (const exception& e) { cerr << "Invalid PT diagnostic interval: " << e.what() << endl; return 1; }
        } else if (strcmp(argv[i], "--load-proposal") == 0 || strcmp(argv[i], "--save-proposal") == 0) {
            const bool load = strcmp(argv[i], "--load-proposal") == 0;
            if (i + 1 >= argc || !*argv[i+1]) { cerr << argv[i] << " requires a filename" << endl; return 1; }
            (load ? proposalInput : proposalOutput) = argv[++i];
        } else if (strcmp(argv[i], "--proposal-warmup") == 0) {
            if (i + 1 >= argc) { cerr << "--proposal-warmup requires a nonnegative integer" << endl; return 1; }
            try {
                string value = argv[++i]; size_t used = 0;
                long long n = stoll(value, &used);
                if (used != value.size() || n < 0 || n > numeric_limits<int>::max())
                    throw invalid_argument("out of range");
                proposalWarmup = static_cast<unsigned>(n);
            } catch (const exception& e) {
                cerr << "Invalid proposal warmup: " << e.what() << endl; return 1;
            }
        } else if (strcmp(argv[i], "--pt-full-priors") == 0) {
            ptFullPriors = true;
        } else if (strcmp(argv[i], "--pt-temperatures") == 0) {
            if (i + 1 >= argc) { cerr << "--pt-temperatures requires a comma-separated list" << endl; return 1; }
            try {
                string list = argv[++i], token;
                if (list.empty() || list.back() == ',') throw invalid_argument("empty temperature");
                istringstream in(list);
                temperatures.clear();
                while (getline(in, token, ',')) {
                    size_t used = 0;
                    double t = stod(token, &used);
                    if (used != token.size() || !isfinite(t) || t < 1. ||
                        (temperatures.empty() ? t != 1. : t <= temperatures.back()))
                        throw invalid_argument("temperatures must increase strictly from 1");
                    temperatures.push_back(t);
                }
                if (temperatures.size() < 2) throw invalid_argument("at least two temperatures required");
            } catch (const exception& e) {
                cerr << "Invalid PT temperatures: " << e.what() << endl; return 1;
            }
        } else if (strcmp(argv[i], "--pt-swap-interval") == 0) {
            if (i + 1 >= argc) { cerr << "--pt-swap-interval requires a positive integer" << endl; return 1; }
            try {
                string value = argv[++i]; size_t used = 0;
                long long n = stoll(value, &used);
                if (used != value.size() || n <= 0 || n > numeric_limits<unsigned>::max())
                    throw invalid_argument("out of range");
                swapInterval = static_cast<unsigned>(n);
            } catch (const exception& e) {
                cerr << "Invalid PT swap interval: " << e.what() << endl; return 1;
            }
        } else if (strcmp(argv[i], "--outfile") == 0 && i + 1 < argc) {
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
        } else if (strcmp(argv[i], "--Univariate") == 0) {
            Univariate = true;
        } else if (strcmp(argv[i], "--indSU3") == 0) {
            flagIndSU3 = true;
        } else if (strcmp(argv[i], "--su3_sigma") == 0 && i + 1 < argc) {
            su3_sigma = atof(argv[++i]);
        } else if (strcmp(argv[i], "--su3_weight") == 0 && i + 1 < argc) {
            su3_weight = argv[++i];
        } else if (strcmp(argv[i], "--gaussianCKM") == 0) {
            flagGaussianCKM = true;
        } else if (strcmp(argv[i], "--positiveG2tP") == 0) {
            flagPositiveG2tP = true;
        } else if (strcmp(argv[i], "--positiveImEA2tP") == 0) {
            flagPositiveImEA2tP = true;
            } else if (strcmp(argv[i], "--smallE2tP") == 0) {
            flagSmallE2tP = true;
        } else if (strcmp(argv[i], "--positiveEA1P") == 0) {
            flagPositiveEA1P = true;
        } else if (strcmp(argv[i], "--smallEA1P") == 0) {
            flagSmallEA1P = true;
        } else if (strcmp(argv[i], "--minimal") == 0) {
            flagMinimal = true;
        } else if (strcmp(argv[i], "--positiveE2tV") == 0) {
            flagPositiveE2tV = true;
        } else if (strcmp(argv[i], "--positiveG2tV") == 0) {
            flagPositiveG2tV = true;
        } else if (strcmp(argv[i], "--positiveG1tD") == 0) {
            flagPositiveG1tD = true;
        } else if (strcmp(argv[i], "--positiveA1D") == 0) {
            flagPositiveA1D = true;
        } else if (strcmp(argv[i], "--help") == 0) {
            cout << "Usage: " << argv[0] << " [options]\n"
                 << "Options:\n"
                 << "  --outfile <filename>            Specify output file name prefix\n"
                 << "  --flagBJPSIV                    Enable BJPSIV flag\n"
                 << "  --flagBJPSIP                    Enable BJPSIP flag\n"
                 << "  --flagBDDb                     Enable BDDb flag\n"
                 << "  --positiveG2tP                  Use positive G2t_scd_BJPSIP values\n"
                 << "  --positiveEA1P                  Use positive EA1_sdcd_BJPSIP values\n"
                 << "  --positiveImEA2tP               Use positive Im(EA2t_ccds_BJPSIP) values\n"
                 << "  --smallE2tP                     Use small E2t_ccdd_BJPSIP values\n"
                 << "  --smallEA1P                     Use small EA1_ccdd_BJPSIP values\n"
                 << "  --minimal                       Use minimal set of observables\n"
                    << "  --positiveE2tV                  Use positive E2t_ccdd_BJPSIV values\n"
                    << "  --positiveG2tV                  Use positive G2t_scd_BJPSIV values\n"
                    << "  --positiveG1tD                  Use positive G1t_scd_BDDb values\n"
                    << "  --positiveA1D                   Use positive A1_sccd_BDDb values\n"
                 << "  --dsu3_limit <value>            Set dsu3 limit (default: 0.2)\n"
                 << "  --ewp_limit <value>             Set ewp limit (default: 0.0)\n"
                 << "  --nPreRun <value>               Set number of pre-run iterations (default: 100000)\n"
                 << "  --nRun <value>                  Set number of run iterations (default: 10000)\n"
                 << "  --nPreRunFactorized <value>     Set number of factorized pre-run iterations (default: 10000)\n"
                  "  --nUpdateMax value            Set max number of update iterations (default: 100)\n"
                  "  --nChains value               Set number of chains (default: 10)\n"
                  "  --save-proposal <file>         Save tuned proposal after pre-run (default: <outfile>proposal.txt)\n"
                  "  --load-proposal <file>         Reuse proposal and replica positions; skip adaptation\n"
                  "  --proposal-warmup <value>      Discard sweeps on load (default: 1000; 0 allowed)\n"
                  "  --pt-temperatures <1,2,...>    Enable PT with explicit temperatures (likelihood only)\n"
                  "  --pt-full-priors              Use unsplit indSU3 prior ranges; override positive/small flags\n"
                  "                                Requires --indSU3 and --pt-temperatures\n"
                  "  --pt-diagnostics-interval N    Save diagnostic blocks every N sweeps (default: 1000)\n"
                  "  --pt-trace-interval N          State trace stride (default: 100; 0 disables traces)\n"
                  "  --no-pt-diagnostics            Disable streaming diagnostics (enabled for PT by default)\n"
                  "  --pt-swap-interval <value>     Sweeps between exchanges (default: 10)\n"
                  "                                With PT, --nChains counts independent ladders\n"
                  "  --Univariate                   Enable univariate mode\n"                  "  --indSU3                        Use goldenmodesB_indSU3 (independent SU(3) breaking)\n"
                  "  --su3_weight <value>            SU(3) weight for indSU3 mode (default: \"abs\")\n"
                  "  --su3_sigma <value>             SU(3) sigma for indSU3 mode (default: free BAT parameter)\n"
                  "  --gaussianCKM                   Use Gaussian priors on CKM parameters (indSU3 mode only)\n"
                  "  --help                          Show this help message\n";
            return 0;
        }
    }

    if (diagnosticOptions && temperatures.empty()) {
        cerr << "PT diagnostic options require --pt-temperatures" << endl; return 1;
    }
    if (ptFullPriors && (!flagIndSU3 || temperatures.empty())) {
        cerr << "--pt-full-priors requires --indSU3 and --pt-temperatures" << endl;
        return 1;
    }
    if (nChains <= 0 || nIterationsRun <= 0 || nIterationsPreRun <= 0 ||
        nIterationsPreRunFactorized < 0 || nIterationsUpdateMax <= 0 ||
        (!temperatures.empty() && temperatures.size() > numeric_limits<unsigned>::max() / unsigned(nChains))) {
        cerr << "Invalid chain/iteration counts" << endl; return 1;
    }
    // Initialize MPI first so rank is known before any output or model construction
    MPI_Init(&argc, &argv);
    int mpi_rank = 0, mpi_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    if (mpi_rank == 0) {
        cout << "Output file name: " << outputFileName << endl;
        cout << "FlagBJPSIV: " << (flagBJPSIV ? "true" : "false") << endl;
        cout << "FlagBJPSIP: " << (flagBJPSIP ? "true" : "false") << endl;
        cout << "FlagBDDb: " << (flagBDDb ? "true" : "false") << endl;
        cout << "dsu3_limit: " << dsu3_limit << endl;
        cout << "ewp_limit: " << ewp_limit << endl;
        cout << "indSU3 mode: " << (flagIndSU3 ? "true" : "false") << endl;
        if (flagIndSU3) {
            cout << "Full prior ranges (PT): " << (ptFullPriors ? "true (positive/small region flags ignored)" : "false") << endl;
            cout << "positive G2t_scd_BJPSIP: " << (flagPositiveG2tP ? "true" : "false") << endl;
            cout << "positive EA1_sdcd_BJPSIP: " << (flagPositiveEA1P ? "true" : "false") << endl;
            cout << "positive Im(EA2t_ccds_BJPSIP): " << (flagPositiveImEA2tP ? "true" : "false") << endl;
            cout << "small E2t_ccdd_BJPSIP: " << (flagSmallE2tP ? "true" : "false") << endl;
            cout << "small EA1_ccdd_BJPSIP: " << (flagSmallEA1P ? "true" : "false") << endl;
            cout << "minimal observables: " << (flagMinimal ? "true" : "false") << endl;
            cout << "su3_sigma: " << (su3_sigma > 0. ? to_string(su3_sigma) : "free") << endl;
            cout << "su3_weight: " << su3_weight << endl;
            cout << "Gaussian CKM prior: " << (flagGaussianCKM ? "true" : "false") << endl;
        }
        cout << "nIterationsPreRun: " << nIterationsPreRun << endl;
        cout << "nIterationsRun: " << nIterationsRun << endl;
        cout << "nIterationsPreRunFactorized: " << nIterationsPreRunFactorized << endl;
        cout << "nIterationsUpdateMax: " << nIterationsUpdateMax << endl;
        cout << "MPI processes: " << mpi_size << endl;
        // Initialize BAT logging (rank 0 only)
        BCLog::OpenLog(outputFileName+"log.txt", BCLog::detail, BCLog::detail);
        BCLog::SetLogLevelScreen(BCLog::summary);
    }

    // Start time measurement
    auto start = chrono::high_resolution_clock::now();

    try {
        if (flagIndSU3) {
            goldenmodesB_indSU3 model(ewp_limit, flagBJPSIP, flagBJPSIV, flagBDDb, su3_sigma, flagGaussianCKM, flagPositiveG2tP, flagPositiveEA1P, flagPositiveImEA2tP, flagSmallE2tP, flagMinimal, flagPositiveE2tV, flagPositiveG2tV, flagSmallEA1P, flagPositiveG1tD, flagPositiveA1D, ptFullPriors);
            model.SetSU3Weight(su3_weight);
            if (mpi_rank != 0) {
                workerLoop(model);
            } else {
                model.procnum = mpi_size;
                cout << "Model (indSU3) constructed correctly, running with "
                     << mpi_size << " MPI process(es)" << endl;
                runModel(model, outputFileName, nChains, nIterationsPreRun, nIterationsRun,
                         nIterationsPreRunFactorized, nIterationsUpdateMax, Univariate, temperatures, swapInterval, proposalInput, proposalOutput, proposalWarmup, ptDiagnostics, diagnosticInterval, traceInterval);
                terminateWorkers(mpi_size, static_cast<int>(model.GetNParameters()) + 1);
            }
        } else {
            goldenmodesB model(dsu3_limit, ewp_limit, flagBJPSIP, flagBJPSIV, flagBDDb);
            if (mpi_rank != 0) {
                workerLoop(model);
            } else {
                model.procnum = mpi_size;
                cout << "Model constructed correctly, running with "
                     << mpi_size << " MPI process(es)" << endl;
                runModel(model, outputFileName, nChains, nIterationsPreRun, nIterationsRun,
                         nIterationsPreRunFactorized, nIterationsUpdateMax, Univariate, temperatures, swapInterval, proposalInput, proposalOutput, proposalWarmup, ptDiagnostics, diagnosticInterval, traceInterval);
                terminateWorkers(mpi_size, static_cast<int>(model.GetNParameters()) + 1);
            }
        }
    } catch (const exception& e) {
        cerr << "Rank " << mpi_rank << ": " << e.what() << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
    if (mpi_rank == 0) {
        BCLog::OutSummary("MCMC analysis completed.");
        BCLog::CloseLog();
        auto endTime = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsedTime = endTime - start;
        cout << "\nExecution time: " << elapsedTime.count() << " seconds" << endl;
    }

    MPI_Finalize();
    return 0;
}
