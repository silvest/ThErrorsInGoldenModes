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

using namespace std;

// Template helper: run MCMC and save results for any model sharing the BCModel+SaveHistograms interface
template <typename Model>
void runModel(Model &model, const string &outputFileName,
              int nChains, int nIterationsPreRun, int nIterationsRun,
              int nIterationsPreRunFactorized, int nIterationsUpdateMax,
              bool Univariate)
{
    model.SetNChains(nChains);
    model.SetNIterationsPreRunMax(nIterationsPreRun);
    model.SetNIterationsRun(nIterationsRun);
    model.SetProposeMultivariate(!Univariate);
    model.SetNIterationsPreRunFactorized(nIterationsPreRunFactorized);
    model.SetNIterationsPreRunCheck(nIterationsUpdateMax);
    model.SetInitialPositionScheme(BCEngineMCMC::kInitRandomPrior);

    for (unsigned i = 0; i < model.GetNParameters(); ++i)
        cout << "Parameter " << i << ": " << model.GetParameter(i).GetName()
             << " [" << model.GetParameter(i).GetLowerLimit()
             << ", " << model.GetParameter(i).GetUpperLimit() << "]" << endl;

    model.MarginalizeAll();
    model.FindModeWithMinuit();
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
    bool flagSU3ReIm = false;
    bool flagGaussianCKM = false;

    double dsu3_limit = 0.2;
    double ewp_limit = 0.0;
    double su3_sigma = -1.0;  // <=0 means free BAT parameter
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
        } else if (strcmp(argv[i], "--Univariate") == 0) {
            Univariate = true;
        } else if (strcmp(argv[i], "--indSU3") == 0) {
            flagIndSU3 = true;
        } else if (strcmp(argv[i], "--su3_sigma") == 0 && i + 1 < argc) {
            su3_sigma = atof(argv[++i]);
        } else if (strcmp(argv[i], "--su3_reIm") == 0) {
            flagSU3ReIm = true;
        } else if (strcmp(argv[i], "--gaussianCKM") == 0) {
            flagGaussianCKM = true;
        } else if (strcmp(argv[i], "--help") == 0) {
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
                  "  --Univariate                   Enable univariate mode\n"                  "  --indSU3                        Use goldenmodesB_indSU3 (independent SU(3) breaking)\n"
                  "  --su3_sigma <value>             SU(3) sigma for indSU3 mode (default: free BAT parameter)\n"
                  "  --su3_reIm                      Use separate Re/Im SU(3) weight instead of |A1-A2|\n"
                  "  --gaussianCKM                   Use Gaussian priors on CKM parameters (indSU3 mode only)\n"
                  "  --help                          Show this help message\n";
            return 0;
        }
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
            cout << "su3_sigma: " << (su3_sigma > 0. ? to_string(su3_sigma) : "free") << endl;
            cout << "su3_reIm weight: " << (flagSU3ReIm ? "true" : "false") << endl;
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

    if (flagIndSU3) {
        goldenmodesB_indSU3 model(ewp_limit, flagBJPSIP, flagBJPSIV, flagBDDb, su3_sigma, flagGaussianCKM);
        model.SetSU3WeightReIm(flagSU3ReIm);
        if (mpi_rank != 0) {
            workerLoop(model);
        } else {
            model.procnum = mpi_size;
            cout << "Model (indSU3) constructed correctly, running with "
                 << mpi_size << " MPI process(es)" << endl;
            runModel(model, outputFileName, nChains, nIterationsPreRun, nIterationsRun,
                     nIterationsPreRunFactorized, nIterationsUpdateMax, Univariate);
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
                     nIterationsPreRunFactorized, nIterationsUpdateMax, Univariate);
            terminateWorkers(mpi_size, static_cast<int>(model.GetNParameters()) + 1);
        }
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
