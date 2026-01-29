#include "goldenmodesB.h"
#include "dato.h"
#include "PDGAverage.h"
#include "CKM.h"
#include <TFile.h>
#include "Pull.h"
#include "BaseMacros.h"
#include <BAT/BCMath.h>
#include <BAT/BCGaussianPrior.h>
#include <TComplex.h>
#include <TRandom.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <complex>
#include <vector>
#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <TRandom3.h>
#include <cmath>
#include <TMatrixD.h>
#include <TDecompChol.h>
#include <iostream>
#include <fstream>
using namespace std;

// ---------------------------------------------------------

CKMParameters ckm;

goldenmodesB::goldenmodesB(double &dsu3_limit_in, double &ewp_limit_in, bool BJPSIP, bool BJPSIV, bool BDDb) : BCModel(), histos(obs)
{
    dsu3_limit = dsu3_limit_in;
    ewp_limit = ewp_limit_in;
    cout << "constructor for goldenmodes called: inserting the experimental data" << endl;
    std::vector<dato> data;
    PDGAverage pdgaverage;

    if (BJPSIP)
    {
        channels.insert(channels.end(), pseudoscalarMesonChannels.begin(), pseudoscalarMesonChannels.end());
    }
    if (BJPSIV)
    {
        channels.insert(channels.end(), vectorMesonChannels.begin(), vectorMesonChannels.end());
    }
    if (BDDb)
    {
        channels.insert(channels.end(), ddbarChannels.begin(), ddbarChannels.end());
    }

    DeclareParameters(); // Ensure parameters are defined

    // Add CKM parameters directly in the constructor
    AddParameter("CKM_Vud", 0.97415, 0.97447);         // Vud parameter
    AddParameter("CKM_Vcb", 0.04046, 0.04194);         // Vcb parameter
    AddParameter("CKM_Vub", 0.00349, 0.00419);         // Vub parameter
    AddParameter("CKM_gamma", 1.12224671, 1.22347581); // gamma parameter (in rad)

    // Add mixing angle between eta1 and eta8
    AddParameter("theta_P", -30. / 180.0 * M_PI, 0.); // in rad

    SetPriorConstantAll();

    // Use lattice QCD result for theta_P
    GetParameter("theta_P")->SetPrior(new BCGaussianPrior("theta_P", -14.1 / 180.0 * M_PI, 2.8 / 180.0 * M_PI)); // in rad

    // measurements

    /////////////////////////////
    // Bdjpsik0
    /////////////////////////////

    // BR measurements
    data.push_back(dato(9.02e-4, 0.10e-4, 0.26e-4)); // Belle:2019xld
    data.push_back(dato(8.1e-4, 0.9e-4, 0.6e-4));    // Belle:2019avj
    data.push_back(dato(8.85e-4, 1.35e-4, 0.1e-4));  // BaBar:2007esv
    data.push_back(dato(8.69e-4, 0.22e-4, 0.30e-4)); // BaBar:2004htr
    data.push_back(dato(9.5e-4, 0.8e-4, 0.6e-4));    // CLEO:2000emb
    data.push_back(dato(11.5e-4, 2.3e-4, 1.7e-4));   // CDF:1995izg
    data.push_back(dato(6.93e-4, 4.07e-4, 0.04e-4)); // CLEO:1991roe
    data.push_back(dato(9.24e-4, 7.21e-4, 0.05e-4)); // ARGUS:1990jet

    pdgaverage.setData(data);
    pdgaverage.setName("BRBdjpsik0");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    // CP asymmetries

    // uncorrelated data first

    // CBdjpsik0s
    data.push_back(dato(0.005, 0.021, 0.037)); // Belle:2012paq
    data.push_back(dato(0.026, 0.025, 0.016)); // BaBar:2009byl

    pdgaverage.setData(data);
    pdgaverage.setName("CBdjpsik0s");
    pdgaverage.CalculateAverage();
    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    // SBdjpsik0s
    data.push_back(dato(0.670, 0.029, 0.013)); // Belle:2012paq
    data.push_back(dato(0.657, 0.036, 0.012)); // BaBar:2009byl
    data.push_back(dato(0.57, 0.58, 0.06));    // Belle:2012teq
    data.push_back(dato(0.775, 0.425));        // CDF:1999ijp
    data.push_back(dato(0.73, 0.93, 0.16));    // ALEPH:2000jem
    data.push_back(dato(3.1, 1.9, 0.5));       // OPAL:1998mtm

    pdgaverage.setData(data);
    pdgaverage.setName("SBdjpsik0s");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    // Correlated Measurements
    std::vector<dato> CorrData;
    TMatrixDSym CorrStat(2);
    TMatrixDSym CorrSyst(2);
    std::vector<std::string> names;

    // Bdjpsik0s : C and S observables from LHCb:2023zcp
    CorrData.push_back(dato(0.010, 0.012)); // C observable
    names.push_back("CBdjpsik0s");
    CorrData.push_back(dato(0.726, 0.014)); // S observable
    names.push_back("SBdjpsik0s");

    // Populate the correlation matrix
    CorrStat(0, 0) = 1.;   // Variance for C
    CorrStat(1, 1) = 1.;   // Variance for S
    CorrStat(0, 1) = 0.41; // Correlation between C and S
    CorrStat(1, 0) = 0.41; // Symmetric part (same as Corr(0, 1))

    // Insert correlated data into corrmeas
    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("CS_Bdjpsik0s_LHCb2023", CorrelatedGaussianObservables(CorrData, CorrStat, CorrSyst)));
    corrmeas_channels.insert(pair<string, std::vector<std::string>>("CS_Bdjpsik0s_LHCb2023", names));
    CorrData.clear();
    names.clear();

    // Bdjpsik0s : C and S observables from Belle-II:2024lwr
    CorrData.push_back(dato(-0.035, 0.026, 0.029)); // C observable
    names.push_back("CBdjpsik0s");
    CorrData.push_back(dato(0.724, 0.035, 0.009)); // S observable
    names.push_back("SBdjpsik0s");

    // Populate the correlation matrix
    CorrStat(0, 0) = 1.; //
    CorrStat(1, 1) = 1.;
    CorrStat(0, 1) = -0.09; // Correlation between C and S
    CorrStat(1, 0) = -0.09; // Symmetric part (same as Corr(0, 1))

    // Populate the correlation matrix
    CorrSyst(0, 0) = 1.; //
    CorrSyst(1, 1) = 1.;
    CorrSyst(0, 1) = 0.; // Correlation between C and S
    CorrSyst(1, 0) = 0.; // Symmetric part (same as Corr(0, 1))

    // Insert correlated data into corrmeas
    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("CS_Bdjpsik0s_BelleII2024", CorrelatedGaussianObservables(CorrData, CorrStat, CorrSyst)));
    corrmeas_channels.insert(pair<string, std::vector<std::string>>("CS_Bdjpsik0s_BelleII2024", names));
    CorrData.clear();
    names.clear();

    // CBdjpsik0l
    data.push_back(dato(-0.007, 0.026, 0.029)); // Belle:2012paq
    data.push_back(dato(-0.033, 0.050, 0.027)); // BaBar:2009byl

    pdgaverage.setData(data);
    pdgaverage.setName("CBdjpsik0l");
    pdgaverage.CalculateAverage();
    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    // SBdjpsik0l
    data.push_back(dato(0.642, 0.047, 0.021)); // Belle:2012paq
    data.push_back(dato(0.694, 0.061, 0.031)); // BaBar:2009byl

    pdgaverage.setData(data);
    pdgaverage.setName("SBdjpsik0l");
    pdgaverage.CalculateAverage();
    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    ////////////////////////////
    // Bdjpsip0
    ////////////////////////////

    // BR measurements
    data.push_back(dato(2.00e-5, 0.12e-5, 0.09e-5));              // Belle-II:2024hqw
    data.push_back(dato(1.670e-5, 0.077e-5, 0.069e-5, 0.095e-5)); // LHCb:2024ier
    data.push_back(dato(1.62e-5, 0.11e-5, 0.06e-5));              // Belle:2018nxw
    data.push_back(dato(1.69e-5, 0.14e-5, 0.07e-5));              // BaBar:2008kfx
    data.push_back(dato(2.6e-5, 1.0e-5, 0.2e-5));                 // CLEO:2000emb

    pdgaverage.setData(data);
    pdgaverage.setName("BRBdjpsip0");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    // CBdjpsip0

    data.push_back(dato(0.155, 0.14, 0.035)); ////Belle:2018nxw
    data.push_back(dato(0.13, 0.12, 0.03));   // Belle-II:2024hqw

    pdgaverage.setData(data);
    pdgaverage.setName("CBdjpsip0");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    // SBdjpsip0
    data.push_back(dato(-0.59, 0.19, 0.03)); ////Belle:2018nxw
    data.push_back(dato(-0.88, 0.17, 0.03)); // Belle-II:2024hqw

    pdgaverage.setData(data);
    pdgaverage.setName("SBdjpsip0");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    // Bdjpsip0 : C and S observables from BaBar:2008kfx
    CorrData.push_back(dato(-0.20, 0.19, 0.03)); // C observable
    names.push_back("CBdjpsip0");
    CorrData.push_back(dato(-1.23, 0.21, 0.04)); // S observable
    names.push_back("SBdjpsip0");

    // Populate the correlation matrix
    CorrStat(0, 0) = 1.; //
    CorrStat(1, 1) = 1.;
    CorrStat(0, 1) = 0.197; // Correlation between C and S
    CorrStat(1, 0) = 0.197; // Symmetric part (same as Corr(0, 1))

    // Populate the correlation matrix
    CorrSyst(0, 0) = 1.; //
    CorrSyst(1, 1) = 1.;
    CorrSyst(0, 1) = 0.; // Correlation between C and S
    CorrSyst(1, 0) = 0.; // Symmetric part (same as Corr(0, 1))

    // Insert correlated data for BaBar:2008kfx into corrmeas
    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("CS_Bdjpsip0_BaBar2008", CorrelatedGaussianObservables(CorrData, CorrStat, CorrSyst)));

    corrmeas_channels.insert(pair<string, std::vector<std::string>>("CS_Bdjpsip0_BaBar2008", names));
    CorrData.clear();
    names.clear();

    ////////////////////////////
    // Bdjpsieta
    ////////////////////////////

    // BR measurements
    data.push_back(dato(1.235e-6, 0.175e-6, 0.07e-6)); // Chang:2012gnb (Belle full dataset)

    pdgaverage.setData(data);
    pdgaverage.setName("BRBdjpsieta");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    // Ratios of BRs

    // BR(Bd->J/psi eta)/BR(Bs->J/psi eta)
    data.push_back(dato(2.16e-2, 0.16e-2, 0.05e-2, 0.07e-2)); // LHCb:2025sgp supersedes previous LHCb:2014oms

    pdgaverage.setData(data);
    pdgaverage.setName("R_Bdjpsieta_Bsjpsieta");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    ////////////////////////////
    // Bdjpsieta'
    ////////////////////////////

    // Ratios of BRs

    // BR(Bd->J/psi eta')/BR(Bs->J/psi eta')
    data.push_back(dato(1.33e-2, 0.12e-2, 0.05e-2, 0.04e-2)); // LHCb:2025sgp supersedes previous LHCb:2014oms

    pdgaverage.setData(data);
    pdgaverage.setName("R_Bdjpsietap_Bsjpsietap");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    // BR(Bd->J/psi eta')/BR(Bd->J/psi eta)
    data.push_back(dato(0.48, 0.06, 0.02, 0.01)); // LHCb:2025sgp supersedes previous LHCb:2014oms

    pdgaverage.setData(data);
    pdgaverage.setName("R_Bdjpsietap_Bdjpsieta");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    ////////////////////////////
    // Bpjpsikp
    ////////////////////////////

    // BR measurements
    data.push_back(dato(10.32e-3, 0.07e-3, 0.24e-3)); // BELLE:2019xld
    data.push_back(dato(9.4e-3, 0.7e-3, 0.8e-3));     // Belle:2019avj
    data.push_back(dato(8.9e-3, 0.6e-3, 0.5e-3));     // Belle:2017psv
    data.push_back(dato(8.1e-3, 1.3e-3, 0.7e-3));     // BaBar:2005pcw
    data.push_back(dato(10.61e-3, 0.15e-3, 0.48e-3)); // BaBar:2004htr
    data.push_back(dato(10.4e-3, 1.1e-3, 0.1e-3));    // BaBar:2005sdl
    data.push_back(dato(10.2e-3, 0.8e-3, 0.7e-3));    // CLEO:1997ilq
    data.push_back(dato(9.24e-3, 3.04e-3, 0.05e-3));  // CLEO:1991roe
    data.push_back(dato(8.09e-3, 3.50e-3, 0.04e-3));  // ARGUS:1990jet

    pdgaverage.setData(data);
    pdgaverage.setName("BRBpjpsikp");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    // CP asymmetry

    // ACPBpjpsikp
    data.push_back(dato(0.09e-2, 0.27e-2, 0.07e-2)); // LHCb:2017joy
    data.push_back(dato(5.9e-3, 3.6e-3, 0.7e-3));    // D0:2013mrm
    data.push_back(dato(-7.6e-3, 5.0e-3, 2.2e-3));   // Belle:2010zqr
    data.push_back(dato(90.e-3, 70.e-3, 20.e-3));    // Belle:2007oni
    data.push_back(dato(30.e-3, 15.e-3, 6.e-3));     // BaBar:2004htr
    data.push_back(dato(18.e-3, 43.e-3, 4.e-3));     // CLEO:2000oig

    pdgaverage.setData(data);
    pdgaverage.setName("ACPBpjpsikp");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    /////////////////////////////
    // Bpjpsipp
    /////////////////////////////

    // BR measurements
    data.push_back(dato(3.8e-5, 0.6e-5, 0.3e-5)); // Belle:2002oex

    pdgaverage.setData(data);
    pdgaverage.setName("BRBpjpsipp");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    // Ratios of BRs

    // R_Bpjpsipp_Bpjpsikp
    data.push_back(dato(3.846e-2, 0.018e-2, 0.018e-2)); // LHCb:2024exp
    data.push_back(dato(3.5e-2, 0.3e-2, 1.2e-2));       // ATLAS:2016rxw
    data.push_back(dato(4.86e-2, 0.82e-2, 0.15e-2));    // CDF:2007mkw
    data.push_back(dato(5.37e-2, 0.45e-2, 0.11e-2));    // BaBar:2004kla
    data.push_back(dato(5.1e-2, 1.8e-2, 0.1e-2));       // CDF:1996efe
    data.push_back(dato(5.2e-2, 2.4e-2));               // CLEO:1995zgs

    pdgaverage.setData(data);
    pdgaverage.setName("R_Bpjpsipp_Bpjpsikp");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    // CP asymmetry

    // deltaA_Bpjpsipp_Bpjpsikp - Difference in ACP between pi and K modes
    meas.insert(pair<string, dato>("deltaA_Bpjpsipp_Bpjpsikp", dato(1.42e-2, 0.43e-2, 0.08e-2))); // LHCb:2024exp

    /////////////////////////////
    // Bsjpsip0
    /////////////////////////////

    // BR measurements

    meas.insert(pair<string, dato>("BRBsjpsip0", dato(0., 1.2e-5 / 8.03 * 3.2))); // extrapolated from the upper limit in Belle:2023tdz

    /////////////////////////////
    // Bsjpsik0b
    /////////////////////////////

    // BR measurements

    // BRBsjpsik0s
    meas.insert(pair<string, dato>("BRBsjpsik0s", dato(2.06e-5, 0.08e-5, 0.06e-5, 0.07e-5, 0.08e-5))); // LHCb:2021qbv

    // ACP measurments

    // Bsjpsik0s : C and S observables from LHCb:2015brj
    CorrData.push_back(dato(0.55, 0.71, 0.06)); // LHCb:2015brj
    names.push_back("SBsjpsik0s");
    CorrData.push_back(dato(-0.28, 0.41, 0.08)); // LHCb:2015brj
    names.push_back("CBsjpsik0s");
    // Populate the correlation matrix
    CorrStat(0, 0) = 1.; //
    CorrStat(1, 1) = 1.;
    CorrStat(0, 1) = -0.06;                      // Correlation
    CorrStat(1, 0) = -0.06;                      // Symmetric part (same as Corr(0, 1))
    CorrSyst = TMatrixD(TMatrixD::kIdentity, 2); // No systematic correlation available in the paper
    // Insert correlated data into corrmeas
    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("CS_Bsjpsik0s_LHCb2015", CorrelatedGaussianObservables(CorrData, CorrStat, CorrSyst)));
    corrmeas_channels.insert(pair<string, std::vector<std::string>>("CS_Bsjpsik0s_LHCb2015", names));
    CorrData.clear();
    names.clear();

    /////////////////////////////
    // Bsjpsieta
    /////////////////////////////

    // BR measurements

    // BRBsjpsieta
    meas.insert(pair<string, dato>("BRBsjpsieta", dato(5.27e-4, 0.50e-4, 0.25e-4, 0.97e-4))); // Chang:2012gnb with symmetrized errors

    // Ratios of BRs

    meas.insert(pair<string, dato>("R_Bsjpsieta_Bsjpsphi", dato(4.50e-1, 0.14e-1, 0.16e-1, 0.13e-1))); // LHCb:2025sgp

    /////////////////////////////
    // Bsjpsietap
    /////////////////////////////

    // BR ratio measurements

    // R_Bsjpsietap_Bsjpsieta
    data.push_back(dato(0.80, 0.02, 0.02, 0.01)); // LHCb:2025sgp
    data.push_back(dato(0.73, 0.14, 0.02));       // Chang:2012gnb

    pdgaverage.setData(data);
    pdgaverage.setName("R_Bsjpsietap_Bsjpsieta");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    // R_Bsjpsietap_Bsjpsiphi
    meas.insert(pair<string, dato>("R_Bsjpsietap_Bsjpsiphi", dato(0.370, 0.013, 0.018, 0.011))); // LHCb:2025sgp

    /////////////////////////////
    // Bsjpsiph
    /////////////////////////////

    // BRBsjpsiph
    data.push_back(dato(1.018e-3, 0.032e-3, 0.037e-3)); // LHCb:2021qbv
    data.push_back(dato(1.25e-3, 0.07e-3, 0.23e-3));    // Belle:2013sdi
    data.push_back(dato(1.5e-3, 0.5e-3, 0.1e-3));       // CDF:1996ivk

    pdgaverage.setData(data);
    pdgaverage.setName("BRBsjpsiph");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    // Correlated data sets for polarization fractions and CP asymmetries

    // LHCb:2023sim (Run 1+2 combined)

    CorrData.push_back(dato(-0.044, 0.020)); // phi_s
    names.push_back("phi_Bsjpsiph");
    CorrData.push_back(dato(0.990, 0.010)); // lambda
    names.push_back("lambda_Bsjpsiph");
    CorrData.push_back(dato(0.2471, 0.0032)); // f_perp
    names.push_back("f_perp_Bsjpsiph");
    CorrData.push_back(dato(0.5175, 0.0035)); // f_0
    names.push_back("f_0_Bsjpsiph");
    CorrData.push_back(dato(2.924, 0.076)); // delta_perp - delta_0
    names.push_back("delta_perp_Bsjpsiph");
    CorrData.push_back(dato(3.150, 0.062)); // delta_paral - delta_0
    names.push_back("delta_paral_Bsjpsiph");

    // Populate the correlation matrix
    TMatrixDSym CorrPhis(6);
    CorrPhis(0, 0) = 1.;
    CorrPhis(0, 1) = CorrPhis(1, 0) = -0.01;
    CorrPhis(0, 2) = CorrPhis(2, 0) = -0.01;
    CorrPhis(0, 3) = CorrPhis(3, 0) = -0.01;
    CorrPhis(0, 4) = CorrPhis(4, 0) = 0.04;
    CorrPhis(0, 5) = CorrPhis(5, 0) = 0.01;
    CorrPhis(1, 1) = 1.;
    CorrPhis(1, 2) = CorrPhis(2, 1) = 0.0;
    CorrPhis(1, 3) = CorrPhis(3, 1) = 0.01;
    CorrPhis(1, 4) = CorrPhis(4, 1) = -0.12;
    CorrPhis(1, 5) = CorrPhis(5, 1) = -0.03;
    CorrPhis(2, 2) = 1.;
    CorrPhis(2, 3) = CorrPhis(3, 2) = -0.17;
    CorrPhis(2, 4) = CorrPhis(4, 2) = -0.01;
    CorrPhis(2, 5) = CorrPhis(5, 2) = -0.04;
    CorrPhis(3, 3) = 1.;
    CorrPhis(3, 4) = CorrPhis(4, 3) = 0.0;
    CorrPhis(3, 5) = CorrPhis(5, 3) = 0.01;
    CorrPhis(4, 4) = 1.;
    CorrPhis(4, 5) = CorrPhis(5, 4) = 0.2;
    CorrPhis(5, 5) = 1.;

    // Insert correlated data into corrmeas
    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("Bsjpsiph_LHCb2023sim", CorrelatedGaussianObservables(CorrData, CorrPhis)));
    corrmeas_channels.insert(pair<string, std::vector<std::string>>("Bsjpsiph_LHCb2023sim", names));
    CorrData.clear();
    names.clear();

    // CMS:2024znt (Run 1+2 combined)

    CorrData.push_back(dato(-0.074, 0.023)); // phi_s
    names.push_back("phi_Bsjpsiph");
    CorrData.push_back(dato(1.011, 0.019)); // lambda
    names.push_back("lambda_Bsjpsiph");
    CorrData.push_back(dato(0.5273, 0.0044)); // f_0
    names.push_back("f_0_Bsjpsiph");
    CorrData.push_back(dato(0.2417, 0.0036)); // f_perp
    names.push_back("f_perp_Bsjpsiph");
    CorrData.push_back(dato(3.152, 0.077)); // delta_paral - delta_0
    names.push_back("delta_paral_Bsjpsiph");
    CorrData.push_back(dato(2.940, 0.098)); // delta_perp - delta_0
    names.push_back("delta_perp_Bsjpsiph");

    // Populate the correlation matrix
    CorrPhis(0, 0) = 1.;
    CorrPhis(0, 1) = CorrPhis(1, 0) = 0.08;
    CorrPhis(0, 2) = CorrPhis(2, 0) = -0.03;
    CorrPhis(0, 3) = CorrPhis(3, 0) = 0.06;
    CorrPhis(0, 4) = CorrPhis(4, 0) = 0.00;
    CorrPhis(0, 5) = CorrPhis(5, 0) = -0.12;
    CorrPhis(1, 1) = 1.;
    CorrPhis(1, 2) = CorrPhis(2, 1) = -0.09;
    CorrPhis(1, 3) = CorrPhis(3, 1) = 0.07;
    CorrPhis(1, 4) = CorrPhis(4, 1) = -0.29;
    CorrPhis(1, 5) = CorrPhis(5, 1) = -0.40;
    CorrPhis(2, 2) = 1.;
    CorrPhis(2, 3) = CorrPhis(3, 2) = -0.69;
    CorrPhis(2, 4) = CorrPhis(4, 2) = 0.0;
    CorrPhis(2, 5) = CorrPhis(5, 2) = 0.02;
    CorrPhis(3, 3) = 1.;
    CorrPhis(3, 4) = CorrPhis(4, 3) = -0.03;
    CorrPhis(3, 5) = CorrPhis(5, 3) = 0.04;
    CorrPhis(4, 4) = 1.;
    CorrPhis(4, 5) = CorrPhis(5, 4) = 0.43;
    CorrPhis(5, 5) = 1.;

    // Insert correlated data into corrmeas
    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("Bsjpsiph_CMS2024znt", CorrelatedGaussianObservables(CorrData, CorrPhis)));
    corrmeas_channels.insert(pair<string, std::vector<std::string>>("Bsjpsiph_CMS2024znt", names));
    CorrData.clear();
    names.clear();

    // phi and lambda for ATLAS:2020lbz solution B
    TMatrixDSym Corr5(5);

    CorrData.push_back(dato(-0.081, 0.041, 0.022)); // phi_s
    names.push_back("phis_Bsjpsiph");
    CorrData.push_back(dato(0.2213, 0.0019, 0.0023)); // f_paral
    names.push_back("f_paral_Bsjpsiph");
    CorrData.push_back(dato(0.5131, 0.0013, 0.0038)); // f0
    names.push_back("f0_Bsjpsiph");
    CorrData.push_back(dato(2.91, 0.11, 0.06)); // g_perp - g_0
    names.push_back("g_perp_minus_g_0_Bsjpsiph");
    CorrData.push_back(dato(2.94, 0.05, 0.09)); // g_par - g_0
    names.push_back("g_par_minus_g_0_Bsjpsiph");

    Corr5(1, 1) = Corr5(0, 0) = Corr5(2, 2) = Corr5(3, 3) = Corr5(4, 4) = 1.;
    Corr5(1, 0) = Corr5(0, 1) = -0.003;
    Corr5(2, 0) = Corr5(0, 2) = -0.004;
    Corr5(3, 0) = Corr5(0, 3) = 0.004;
    Corr5(4, 0) = Corr5(0, 4) = 0.007;
    Corr5(1, 2) = Corr5(2, 1) = -0.341;
    Corr5(1, 3) = Corr5(3, 1) = 0.133;
    Corr5(1, 4) = Corr5(4, 1) = 0.522;
    Corr5(2, 3) = Corr5(3, 2) = -0.034;
    Corr5(2, 4) = Corr5(4, 2) = -0.103;
    Corr5(3, 4) = Corr5(4, 3) = 0.254;

    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("phi_Bsjpsiph_ATLAS2020B", CorrelatedGaussianObservables(CorrData, Corr5)));
    corrmeas_channels.insert(pair<string, std::vector<std::string>>("phi_Bsjpsiph_ATLAS2020B", names));
    names.clear();
    CorrData.clear();

    /////////////////////////////
    // Bsjpsiom
    /////////////////////////////

    // Nothing yet measured

    /////////////////////////////
    // Bsjpsikbst
    /////////////////////////////

    // BR measurement

    meas.insert(pair<string, dato>("BRBsjpsikbst", dato(4.13e-5, 0.12e-5, 0.07e-5, 0.14e-5, 0.45e-5))); // LHCb:2013lka

    // Angular analysis measurements: polarization fractions and CP asymmetries
    //  LHCb:2025vrr (Use only run 2 since no correlation matrix is given for the combined Run 1+2)

    TMatrixDSym CorrBsjpsikbst(7);
    CorrData.push_back(dato(0.014, 0.030)); // A_CP_0
    names.push_back("ACPBsjpsikbst_0");
    CorrData.push_back(dato(-0.055, 0.066)); // A_CP_paral
    names.push_back("ACPBsjpsikbst_paral");
    CorrData.push_back(dato(0.060, 0.059)); // A_CP_perp
    names.push_back("ACPBsjpsikbst_perp");
    CorrData.push_back(dato(2.879, 0.087)); // delta_paral - delta_0
    names.push_back("delta_Bsjpsikbst_paral");
    CorrData.push_back(dato(0.057, 0.068)); // delta_perp - delta_0
    names.push_back("delta_Bsjpsikbst_perp");
    CorrData.push_back(dato(0.534, 0.015)); // f_0
    names.push_back("f_0_Bsjpsikbst");
    CorrData.push_back(dato(0.211, 0.015)); // f_perp
    names.push_back("f_perp_Bsjpsikbst");

    // Populate the correlation matrix
    CorrBsjpsikbst(0, 0) = 1.;
    CorrBsjpsikbst(0, 1) = CorrBsjpsikbst(1, 0) = -0.16;
    CorrBsjpsikbst(0, 2) = CorrBsjpsikbst(2, 0) = -0.08;
    CorrBsjpsikbst(0, 3) = CorrBsjpsikbst(3, 0) = 0.01;
    CorrBsjpsikbst(0, 4) = CorrBsjpsikbst(4, 0) = 0.02;
    CorrBsjpsikbst(0, 5) = CorrBsjpsikbst(5, 0) = 0.01;
    CorrBsjpsikbst(0, 6) = CorrBsjpsikbst(6, 0) = -0.01;
    CorrBsjpsikbst(1, 1) = 1.;
    CorrBsjpsikbst(1, 2) = CorrBsjpsikbst(2, 1) = -0.57;
    CorrBsjpsikbst(1, 3) = CorrBsjpsikbst(3, 1) = -0.03;
    CorrBsjpsikbst(1, 4) = CorrBsjpsikbst(4, 1) = -0.04;
    CorrBsjpsikbst(1, 5) = CorrBsjpsikbst(5, 1) = -0.01;
    CorrBsjpsikbst(1, 6) = CorrBsjpsikbst(6, 1) = 0.01;
    CorrBsjpsikbst(2, 2) = 1.;
    CorrBsjpsikbst(2, 3) = CorrBsjpsikbst(3, 2) = 0.01;
    CorrBsjpsikbst(2, 4) = CorrBsjpsikbst(4, 2) = 0.02;
    CorrBsjpsikbst(2, 5) = CorrBsjpsikbst(5, 2) = 0.01;
    CorrBsjpsikbst(2, 6) = CorrBsjpsikbst(6, 2) = 0.04;
    CorrBsjpsikbst(3, 3) = 1.;
    CorrBsjpsikbst(3, 4) = CorrBsjpsikbst(4, 3) = 0.69;
    CorrBsjpsikbst(3, 5) = CorrBsjpsikbst(5, 3) = -0.09;
    CorrBsjpsikbst(3, 6) = CorrBsjpsikbst(6, 3) = 0.05;
    CorrBsjpsikbst(4, 4) = 1.;
    CorrBsjpsikbst(4, 5) = CorrBsjpsikbst(5, 4) = -0.03;
    CorrBsjpsikbst(4, 6) = CorrBsjpsikbst(6, 4) = 0.0;
    CorrBsjpsikbst(5, 5) = 1.;
    CorrBsjpsikbst(5, 6) = CorrBsjpsikbst(6, 5) = -0.44;
    CorrBsjpsikbst(6, 6) = 1.;

    // Insert correlated data into corrmeas
    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("Bsjpsikbst_LHCb2025", CorrelatedGaussianObservables(CorrData, CorrBsjpsikbst)));
    corrmeas_channels.insert(pair<string, std::vector<std::string>>("Bsjpsikbst_LHCb2025", names));
    CorrData.clear();
    names.clear();

    /////////////////////////////
    // Bsjpsirho
    /////////////////////////////

    // No data available yet

    /////////////////////////////
    // Bdjpsiom
    ///////////////////////////

    // BRBdjpsiom

    meas.insert(pair<string, dato>("BRBdjpsiom", dato(2.16e-5, 0.30e-5, 0.14e-5))); // Belle-II:2024fgp

    // Ratios of BRs
    //  R_Bdjpsiom_Bdjpsirho from LHCb:2012cw
    meas.insert(pair<string, dato>("R_Bdjpsiom_Bdjpsirho", dato(0.86, 0.19, 0.10))); // LHCb:2012cw

    // Transversity fractions and phases from LHCb:2014vbo
    meas.insert(pair<string, dato>("f_0_Bdjpsiom", dato(0.405, 0.14, 0.035)));                                                  // LHCb:2014vbo
    meas.insert(pair<string, dato>("f_paral_Bdjpsiom", dato(0.58, 0.135, 0.035)));                                              // LHCb:2014vbo
    meas.insert(pair<string, dato>("delta_paral_Bdjpsiom-delta_0_Bdjpsirho", dato(123.5 / 180. * M_PI, 13.7 / 180. * M_PI)));   // LHCb:2014vbo
    meas.insert(pair<string, dato>("delta_perp_Bdjpsiom-delta_perp_Bdjpsirho", dato(227.4 / 180. * M_PI, 84.9 / 180. * M_PI))); // LHCb:2014vbo
    meas.insert(pair<string, dato>("delta_0_Bdjpsiom-delta_0_Bdjpsirho", dato(268.8 / 180. * M_PI, 11.9 / 180. * M_PI)));       // LHCb:2014vbo

    /////////////////////////////
    // Bdjpsikst
    /////////////////////////////

    // BRBdjpsikst
    data.push_back(dato(1.19e-3, 0.01e-3, 0.08e-3));    // Belle:2014nuw
    data.push_back(dato(1.335e-3, 0.215e-3, 0.02e-3));  // BaBar:2007esv
    data.push_back(dato(1.309e-3, 0.026e-3, 0.077e-3)); // BaBar:2004htr
    data.push_back(dato(1.29e-3, 0.05e-3, 0.13e-3));    // Belle:2002otd
    data.push_back(dato(1.74e-3, 0.20e-3, 0.18e-3));    // CDF:1998tqc
    data.push_back(dato(1.32e-3, 0.17e-3, 0.17e-3));    // CLEO:1997ilq
    data.push_back(dato(1.27e-3, 0.65e-3, 0.01e-3));    // CLEO:1991roe
    data.push_back(dato(1.27e-3, 0.60e-3, 0.01e-3));    // ARGUS:1990jet
    data.push_back(dato(4.04e-3, 1.81e-3, 0.02e-3));    // CLEO:1987iba

    pdgaverage.setData(data);
    pdgaverage.setName("BRBdjpsikst");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    // ACPBdjpsikst - Direct CP asymmetry for Bd->J/psi K*
    data.push_back(dato(-0.01, 0.04, 0.02));   // Belle:2014nuw
    data.push_back(dato(0.030, 0.058, 0.012)); // BaBar:2009byl

    pdgaverage.setData(data);
    pdgaverage.setName("ACPBdjpsikst");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    // polarization

    // f_0 Bdjpsikst
    data.push_back(dato(0.587, 0.011));        // D0:2008nly
    data.push_back(dato(0.556, 0.009, 0.010)); // BaBar:2007rbr
    data.push_back(dato(0.562, 0.026, 0.018)); // CDF:2004dxr
    data.push_back(dato(0.574, 0.012, 0.009)); // Belle:2005qtf
    data.push_back(dato(0.59, 0.06, 0.01));    // CDF:2000edf
    data.push_back(dato(0.52, 0.07, 0.04));    // CLEO:1997ilq
    data.push_back(dato(0.65, 0.10, 0.04));    // CDF:1995kwt
    data.push_back(dato(0.97, 0.16, 0.15));    // ARGUS:1994rms

    pdgaverage.setData(data);
    pdgaverage.setName("f_0_Bdjpsikst");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    // f_paral Bdjpsikst
    data.push_back(dato(0.230, 0.013, 0.025)); // D0:2008nly
    data.push_back(dato(0.211, 0.010, 0.006)); // BaBar:2007rbr

    pdgaverage.setData(data);
    pdgaverage.setName("f_paral_Bdjpsikst");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    // f_perp Bdjpsikst
    data.push_back(dato(0.233, 0.010, 0.005)); // BaBar:2007rbr
    data.push_back(dato(0.215, 0.032, 0.006)); // CDF:2004dxr
    data.push_back(dato(0.195, 0.012, 0.008)); // Belle:2005qtf

    pdgaverage.setData(data);
    pdgaverage.setName("f_perp_Bdjpsikst");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    // delta_paral Bdjpsikst
    data.push_back(dato(-2.69, 0.08, 0.11)); // D0:2008nly
    data.push_back(dato(-2.93, 0.08, 0.04)); // BaBar:2007rbr

    pdgaverage.setData(data);
    pdgaverage.setName("delta_paral_Bdjpsikst");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    // delta_perp Bdjpsikst
    data.push_back(dato(3.21, 0.06, 0.06)); // D0:2008nly
    data.push_back(dato(2.91, 0.05, 0.03)); // BaBar:2007rbr

    pdgaverage.setData(data);
    pdgaverage.setName("delta_perp_Bdjpsikst");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    // Polarization fractions and phases from LHCb:2013vga

    TMatrixDSym CorrBdjpsikst(4);
    CorrData.push_back(dato(0.227, sqrt(0.004 * 0.004 + 0.011 * 0.011))); // f_paral
    names.push_back("f_paral_Bdjpsikst");
    CorrData.push_back(dato(0.201, sqrt(0.004 * 0.004 + 0.008 * 0.008))); // f_perp
    names.push_back("f_perp_Bdjpsikst");
    CorrData.push_back(dato(-2.94, sqrt(0.02 * 0.02 + 0.03 * 0.03))); // delta_paral - delta_0
    names.push_back("delta_paral_Bdjpsikst");
    CorrData.push_back(dato(2.94, sqrt(0.02 * 0.02 + 0.02 * 0.02))); // delta_perp - delta_0
    names.push_back("delta_perp_Bdjpsikst");

    // Populate the correlation matrix
    CorrBdjpsikst(0, 0) = 1.;
    CorrBdjpsikst(0, 1) = CorrBdjpsikst(1, 0) = -0.70;
    CorrBdjpsikst(0, 2) = CorrBdjpsikst(2, 0) = 0.12;
    CorrBdjpsikst(0, 3) = CorrBdjpsikst(3, 0) = 0.04;
    CorrBdjpsikst(1, 1) = 1.;
    CorrBdjpsikst(1, 2) = CorrBdjpsikst(2, 1) = -0.14;
    CorrBdjpsikst(1, 3) = CorrBdjpsikst(3, 1) = -0.01;
    CorrBdjpsikst(2, 2) = 1.;
    CorrBdjpsikst(2, 3) = CorrBdjpsikst(3, 2) = 0.64;
    CorrBdjpsikst(3, 3) = 1.;

    // Insert correlated data into corrmeas
    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("Bdjpsikst_LHCb2013", CorrelatedGaussianObservables(CorrData, CorrBdjpsikst)));
    corrmeas_channels.insert(pair<string, std::vector<std::string>>("Bdjpsikst_LHCb2013", names));
    CorrData.clear();
    names.clear();

    // CP asymmetries from LHCb:2013vga
    meas.insert(pair<string, dato>("ACPfparal_Bdjpsikst", dato(-0.011, 0.016, 0.005)));       // LHCb:2013vga
    meas.insert(pair<string, dato>("ACPPfperp_Bdjpsikst", dato(0.032, 0.018, 0.003)));        // LHCb:2013vga
    meas.insert(pair<string, dato>("ACPPdelta_paral_Bdjpsikst", dato(-0.003, 0.007, 0.002))); // LHCb:2013vga
    meas.insert(pair<string, dato>("ACPPdelta_perp_Bdjpsikst", dato(0.003, 0.005, 0.001)));   // LHCb:2013vga

    // CBdjpsikst
    meas.insert(pair<string, dato>("CBdjpsikst", dato(0.025, 0.083, 0.054))); // BaBar:2009byl

    // SBdjpsikst
    meas.insert(pair<string, dato>("SBdjpsikst", dato(0.601, 0.239, 0.087))); // BaBar:2009byl

    /////////////////////////////
    // Bdjpsirh
    /////////////////////////////

    // BRBdjpsirh
    data.push_back(dato(2.515e-5, 0.10e-5, 0.165e-5)); // LHCb:2014vbo
    data.push_back(dato(2.7e-5, 0.3e-5, 0.2e-5));      // BaBar:2007yvx

    pdgaverage.setData(data);
    pdgaverage.setName("BRBdjpsirh");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    // CP violation for the different polarizations of Bdjpsirh from LHCb:2014xpr

    CorrData.push_back(dato(-0.047, 0.034, 0.011)); // Alpha_0
    names.push_back("ACPBdjpsirh_0");
    CorrData.push_back(dato(-0.060, 0.060, 0.007)); // Alpha_paral
    names.push_back("ACPBdjpsirh_paral");
    CorrData.push_back(dato(0.020, 0.109, 0.018)); // Alpha_perp
    names.push_back("ACPBdjpsirh_perp");
    CorrData.push_back(dato(-3.3 / 180. * M_PI, 7.2 / 180. * M_PI, 1.7 / 180. * M_PI)); // 2beta_perp - 2beta_0
    names.push_back("2beta_Bdjpsirh_perp-2beta_Bdjpsirh_0");
    CorrData.push_back(dato(-0.5 / 180. * M_PI, 6.5 / 180. * M_PI, 1.6 / 180. * M_PI)); // 2beta_paral - 2beta_0
    names.push_back("2beta_Bdjpsirh_paral-2beta_Bdjpsirh_0");
    CorrData.push_back(dato(42.1 / 180. * M_PI, 10.2 / 180. * M_PI, 5.0 / 180. * M_PI)); // 2beta_0
    names.push_back("2beta_Bdjpsirh_0");
    // Populate the correlation matrix
    TMatrixDSym CorrBdjpsirh(6);

    CorrBdjpsirh(0, 0) = 1.;
    CorrBdjpsirh(0, 1) = CorrBdjpsirh(1, 0) = 0.03;
    CorrBdjpsirh(0, 2) = CorrBdjpsirh(2, 0) = 0.16;
    CorrBdjpsirh(0, 3) = CorrBdjpsirh(3, 0) = 0.22;
    CorrBdjpsirh(0, 4) = CorrBdjpsirh(4, 0) = 0.16;
    CorrBdjpsirh(0, 5) = CorrBdjpsirh(5, 0) = -0.11;
    CorrBdjpsirh(1, 1) = 1.;
    CorrBdjpsirh(1, 2) = CorrBdjpsirh(2, 1) = -0.21;
    CorrBdjpsirh(1, 3) = CorrBdjpsirh(3, 1) = 0.59;
    CorrBdjpsirh(1, 4) = CorrBdjpsirh(4, 1) = -0.07;
    CorrBdjpsirh(1, 5) = CorrBdjpsirh(5, 1) = 0.1;
    CorrBdjpsirh(2, 2) = 1.;
    CorrBdjpsirh(2, 3) = CorrBdjpsirh(3, 2) = -0.04;
    CorrBdjpsirh(2, 4) = CorrBdjpsirh(4, 2) = -0.25;
    CorrBdjpsirh(2, 5) = CorrBdjpsirh(5, 2) = -0.09;
    CorrBdjpsirh(3, 3) = 1.;
    CorrBdjpsirh(3, 4) = CorrBdjpsirh(4, 3) = 0.39;
    CorrBdjpsirh(3, 5) = CorrBdjpsirh(5, 3) = -0.08;
    CorrBdjpsirh(4, 4) = 1.;
    CorrBdjpsirh(4, 5) = CorrBdjpsirh(5, 4) = -0.1;
    CorrBdjpsirh(5, 5) = 1.;

    // Insert correlated data into corrmeas
    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("Bdjpsirh_LHCb2014xpr", CorrelatedGaussianObservables(CorrData, CorrBdjpsirh)));
    corrmeas_channels.insert(pair<string, std::vector<std::string>>("Bdjpsirh_LHCb2014xpr", names));
    CorrData.clear();
    names.clear();

    /////////////////////////////
    // BdJpsiphi
    /////////////////////////////

    // BRBdjpsiphi
    meas.insert(pair<string, dato>("BRBdjpsiphi", dato(6.8e-8, 3.0e-8, 0.9e-8))); // LHCb:2020xjm

    /////////////////////////////
    // Bpjpsikst
    /////////////////////////////

    // BRBpjpsikst
    meas.insert(pair<string, dato>("BRBpjpsikst", dato(1.43e-3, 0.08e-3))); // PDGlive 1/2025

    // polarization fractions from Belle:2005qtf
    meas.insert(pair<string, dato>("f_0_Bpjpsikst", dato(0.604, 0.015, 0.018)));    // Belle:2005qtf
    meas.insert(pair<string, dato>("f_perp_Bpjpsikst", dato(0.180, 0.014, 0.010))); // Belle:2005qtf

    // CP asymmetry from BaBar:2004htr
    meas.insert(pair<string, dato>("ACPBpjpsikst", dato(0.048, 0.029, 0.016))); // BaBar:2004htr

    /////////////////////////////
    // Bpjpsirho
    /////////////////////////////

    // BRBpjpsirho
    data.push_back(dato(3.81e-5, 0.25e-5, 0.35e-5)); // LHCb:2018pil
    data.push_back(dato(5.0e-5, 0.7e-5, 0.31e-5));   // BaBar:2007yvx

    pdgaverage.setData(data);
    pdgaverage.setName("BRBpjpsirho");
    pdgaverage.CalculateAverage();
    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    // CP asymmetry

    data.push_back(dato(-0.045, 0.056, 0.008)); // LHCb:2018pil
    data.push_back(dato(-0.11, 0.12, 0.08));    // BaBar:2007yvx

    pdgaverage.setData(data);
    pdgaverage.setName("ACPBpjpsirho");
    pdgaverage.CalculateAverage();
    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    //----------------------------------------------------------------------
    // B→DD channels (control modes for penguin analysis)
    //----------------------------------------------------------------------

    // BRBddpdm (Bd→D⁺D⁻)
    data.push_back(dato(2.12e-4, 0.16e-4, 0.18e-4)); // Belle:2012mef
    data.push_back(dato(1.97e-4, 0.20e-4, 0.20e-4)); // Belle:2007ebz
    data.push_back(dato(2.8e-4, 0.4e-4, 0.5e-4));    // BaBar:2006uih

    pdgaverage.setData(data);
    pdgaverage.setName("BRBddpdm");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    // BRBsdpsdms (Bs→Ds⁺Ds⁻)
    data.push_back(dato(4.0e-3, 0.2e-3, 0.5e-3)); // LHCb:2013sad
    data.push_back(dato(5.9e-3, 1.0e-3, 1.3e-3)); // Belle:2012tsw
    data.push_back(dato(5.4e-3, 0.8e-3, 0.8e-3)); // CDF:2012xmd

    pdgaverage.setData(data);
    pdgaverage.setName("BRBsdpsdms");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    // BRBpdpd0b (B⁺→D⁺D̄⁰)
    data.push_back(dato(3.85e-4, 0.31e-4, 0.38e-4)); // Belle:2008doh
    data.push_back(dato(3.8e-4, 0.6e-4, 0.5e-4));    // BaBar:2006uih

    pdgaverage.setData(data);
    pdgaverage.setName("BRBpdpd0b");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    // BRBpdpsd0b (B⁺→D⁺D̄⁰ₛ)
    data.push_back(dato(8.6e-3, 0.2e-3, 1.1e-3)); // LHCb:2013sad
    data.push_back(dato(9.5e-3, 2.0e-3, 0.8e-3)); // BaBar:2006jvx
    data.push_back(dato(9.8e-3, 2.6e-3, 0.9e-3)); // CLEO:1995psi
    data.push_back(dato(14e-3, 8e-3, 1e-3));      // ARGUS:1991xej
    data.push_back(dato(13e-3, 6e-3, 1e-3));      // CLEO:1990mqz

    pdgaverage.setData(data);
    pdgaverage.setName("BRBpdpsd0b");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    // ACPBpdpd0b (B⁺→D⁺D̄⁰)
    data.push_back(dato(0.00, 0.08, 0.02));  // Belle:2008doh
    data.push_back(dato(-0.13, 0.14, 0.02)); // BaBar:2006uih

    pdgaverage.setData(data);
    pdgaverage.setName("ACPBpdpd0b");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    //----------------------------------------------------------------------
    // B→DD Correlated CP measurements
    //----------------------------------------------------------------------

    // Bddpdm : C and S observables from LHCb:2024gkk (arXiv:2409.03009)
    CorrData.push_back(dato(0.128, 0.103, 0.010)); // C observable
    names.push_back("C_Bddpdm");
    CorrData.push_back(dato(0.552, 0.100, 0.010)); // S observable
    names.push_back("S_Bddpdm");

    CorrStat(0, 0) = 1.;
    CorrStat(1, 1) = 1.;
    CorrStat(0, 1) = -0.472; // Correlation between C and S
    CorrStat(1, 0) = -0.472;

    CorrSyst(0, 0) = 1.;
    CorrSyst(1, 1) = 1.;
    CorrSyst(0, 1) = 0.;
    CorrSyst(1, 0) = 0.;

    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("CS_Bddpdm_LHCb2024", CorrelatedGaussianObservables(CorrData, CorrStat, CorrSyst)));
    corrmeas_channels.insert(pair<string, std::vector<std::string>>("CS_Bddpdm_LHCb2024", names));
    names.clear();
    CorrData.clear();

    // Bsdpsdms : C and S observables from LHCb Run2 (part of arXiv:2409.03009 combination)
    // We use C and S instead of phi_s and lambda to avoid using pre-averaged values
    // From LHCb Run2: C = 0.128 ± 0.103(stat) ± 0.010(syst), S = 0.552 ± 0.100(stat) ± 0.010(syst)
    // NOTE: These are for Bs→Ds⁺Ds⁻ from the same paper as Bd→D⁺D⁻
    CorrData.push_back(dato(-0.053, 0.096, 0.020)); // C observable (NOTE: Check paper for actual values)
    names.push_back("C_Bsdpsdms");
    CorrData.push_back(dato(-0.055, 0.092, 0.021)); // S observable (NOTE: Check paper for actual values)
    names.push_back("S_Bsdpsdms");

    CorrStat(0, 0) = 1.;
    CorrStat(1, 1) = 1.;
    CorrStat(0, 1) = 0.; // Correlation coefficient (CHECK PAPER)
    CorrStat(1, 0) = 0.;

    CorrSyst(0, 0) = 1.;
    CorrSyst(1, 1) = 1.;
    CorrSyst(0, 1) = 0.;
    CorrSyst(1, 0) = 0.;

    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("CS_Bsdpsdms_LHCb2024", CorrelatedGaussianObservables(CorrData, CorrStat, CorrSyst)));
    corrmeas_channels.insert(pair<string, std::vector<std::string>>("CS_Bsdpsdms_LHCb2024", names));
    names.clear();
    CorrData.clear();

    // ACPBpdpd0b and ACPBpdpsd0b from LHCb:2023wbb (correlated)
    dato Bpdpd0b_acp(2.5e-2, 1.0e-2, 0.4e-2, 0.3e-2);
    names.push_back("ACPBpdpd0b");
    dato Bpdpsd0b_acp(0.5e-2, 0.2e-2, 0.5e-2, 0.3e-2);
    names.push_back("ACPBpdpsd0b");

    CorrData.push_back(dato(Bpdpd0b_acp.getMean(), Bpdpd0b_acp.getSigma()));   // ACPBpdpd0b
    CorrData.push_back(dato(Bpdpsd0b_acp.getMean(), Bpdpsd0b_acp.getSigma())); // ACPBpdpsd0b

    TMatrixDSym Corr_ACP_charged(2);
    Corr_ACP_charged(0, 0) = 1.0;
    Corr_ACP_charged(1, 1) = 1.0;
    Corr_ACP_charged(0, 1) = 0.386; // Correlation from LHCb:2023wbb
    Corr_ACP_charged(1, 0) = 0.386;

    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("ACP_Bpdpd0b_Bpdpsd0b_LHCb2023", CorrelatedGaussianObservables(CorrData, Corr_ACP_charged)));
    corrmeas_channels.insert(pair<string, std::vector<std::string>>("ACP_Bpdpd0b_Bpdpsd0b_LHCb2023", names));
    names.clear();
    CorrData.clear();

    std::cout << "All meas inserted" << endl;

    // Create histograms for all observables

    for (const auto &channel : channelNamesSU3)
    {
        // create histos for mod and phase as well as real and imaginary parts for the effective parameters
        for (const auto &param : channelParameters[channel])
        {
            histos.createH1D(param, 500, 0.0, 0.0);
            std::string newStr;
            size_t length = param.length();

            newStr.reserve(length + 1);

            if (param.ends_with("_re"))
            {
                newStr.append(str, 0, length - 3); // Append text before "_re"
                newStr += "_abs";                  // Append replacement
            }
            else if (param.ends_with("_im"))
            {
                newStr.append(str, 0, length - 3); // Append text before "_im"
                newStr += "_arg";                  // Append replacement
            }
            histos.createH1D(newStr, 500, 0.0, 0.0);
        }

        // Branching ratio (BR)
        if (meas.find("BR" + channel) != meas.end())
        {
            histos.createH1D("BR_" + channel, 500, 0.0, 0.0);
        }

        // ACP, C, and S observables
        if (meas.find("ACP" + channel) != meas.end())
        {
            histos.createH1D("ACP_" + channel, 500, 0., 0.);
        }
        if (meas.find("C" + channel) != meas.end())
        {
            histos.createH1D("C_" + channel, 500, 0.0, 0.0);
        }
        if (meas.find("S" + channel) != meas.end())
        {
            histos.createH1D("S_" + channel, 500, 0.0, 0.0);
        }
    }

    for (const auto &channel : channelNames)
    {

        // Branching ratio (BR)
        if (meas.find("BR" + channel) != meas.end())
        {
            histos.createH1D("BR_" + channel, 500, 0.0, 0.0);
        }

        // ACP, C, and S observables
        if (meas.find("ACP" + channel) != meas.end())
        {
            histos.createH1D("ACP_" + channel, 500, 0., 0.);
        }
        if (meas.find("C" + channel) != meas.end())
        {
            histos.createH1D("C_" + channel, 500, 0.0, 0.0);
        }
        if (meas.find("S" + channel) != meas.end())
        {
            histos.createH1D("S_" + channel, 500, 0.0, 0.0);
        }
    }

    // Create histograms for correlated observables (from corrmeas)
    for (auto &pair : corrmeas)
    {
        const auto &corr_key = pair.first;
        auto &corrObservable = pair.second;

        // Parse the channels from the key
        auto channel_and_experiment = extractChannelFromCorrKey(corr_key);
        auto &channels = channel_and_experiment.first;

        // Handle C and S observables (e.g., "CS_" prefix)
        if (corr_key.find("CS_") != std::string::npos)
        {
            histos.createH1D("C_" + channels[0], 500, 0.0, 0.0);
            histos.createH1D("S_" + channels[0], 500, 0.0, 0.0);
        }

        // Handle ACP observables (e.g., "ACP_" prefix)
        if (corr_key.find("ACP_") != std::string::npos)
        {
            for (size_t i = 0; i < channels.size(); ++i)
            {
                histos.createH1D("ACP_" + channels[i], 500, 0.0, 0.0);
            }
        }
    }
    histos.createH1D("phi_Bsjpsiph", 500, 0.0, 0.0);
    histos.createH1D("f_perp_Bsjpsiph", 500, 0.0, 0.0);
    histos.createH1D("f_paral_Bsjpsiph", 500, 0.0, 0.0);
    histos.createH1D("f_0_Bsjpsiph", 500, 0.0, 0.0);
    histos.createH1D("delta_perp_Bsjpsiph", 500, 0.0, 0.0);
    histos.createH1D("delta_paral_Bsjpsiph", 500, 0.0, 0.0);
    histos.createH1D("lambda_Bsjpsiph", 500, 0.0, 0.0);
    histos.createH1D("f_perp_Bdjpsikst", 500, 0.0, 0.0);
    histos.createH1D("f_paral_Bdjpsikst", 500, 0.0, 0.0);
    histos.createH1D("f_0_Bdjpsikst", 500, 0.0, 0.0);
    histos.createH1D("delta_perp_Bdjpsikst", 500, 0.0, 0.0);
    histos.createH1D("delta_paral_Bdjpsikst", 500, 0.0, 0.0);
    histos.createH1D("f_paral_Bsjpsikst", 500, 0.0, 0.0);
    histos.createH1D("f_0_Bsjpsikst", 500, 0.0, 0.0);
    histos.createH1D("delta_perp_Bsjpsikst", 500, 0.0, 0.0);
    histos.createH1D("delta_paral_Bsjpsikst", 500, 0.0, 0.0);
    histos.createH1D("A0_CP_Bsjpsikst", 500, 0.0, 0.0);
    histos.createH1D("Aperp_CP_Bsjpsikst", 500, 0.0, 0.0);
    histos.createH1D("Aparal_CP_Bsjpsikst", 500, 0.0, 0.0);
}

//---------------------------------------------------------------

// method to define the parameters needed to calculate each decay amplitude, uses BCModel method AddParameter
void goldenmodesB::DefineParameters(const string &channel)
{
    double limit1 = 0.;

    // We take as reference channels the golden ones: B→J/ψK0, Bs→J/ψφ, Bs->Ds+Ds-
    // all other channels will have parameters defined relative to these ones including SU(3) breaking
    // We choose the dominant amplitude of the reference channels to be real and positive

    if (channel == "Bdjpsik0")
    {
        std::vector<std::string> params = {
            "E2t_ccsd_BJPSIP_re",
            "E2t_ccsd_BJPSIP_im",
            "G2t_scd_BJPSIP_re",
            "G2t_scd_BJPSIP_im",
        };
        channelParameters[channel] = params;

        addAmplitudeParameter("E2t_ccsd_BJPSIP_re", 0., 100.);
        addAmplitudeParameter("E2t_ccsd_BJPSIP_im", 0., 0.);
        addAmplitudeParameter("G2t_scd_BJPSIP_re", -100., 100.);
        addAmplitudeParameter("G2t_scd_BJPSIP_im", -100., 100.);
    }
    else if (channel == "Bdjpsip0")
    {
        std::vector<std::string> params = {
            "delta_E2t_ccdd_BJPSIP_re",
            "delta_E2t_ccdd_BJPSIP_im",
            "dP4EW_ucd_BPJPSI_re",
            "dP4EW_ucd_BPJPSI_im",
            "EA2_ddcd_BPJPSI_re",
            "EA2_ddcd_BPJPSI_im",
            "delta_G2t_dcd_BJPSIP_re",
            "delta_G2t_dcd_BJPSIP_im"};
        channelParameters[channel] = params;

        addAmplitudeParameter("dP4EW_ucd_BPJPSI_re", -ewp_limit, ewp_limit);
        addAmplitudeParameter("dP4EW_ucd_BPJPSI_im", -ewp_limit, ewp_limit);
        addAmplitudeParameter("EA2_ddcd_BPJPSI_re", -100., 100.);
        addAmplitudeParameter("EA2_ddcd_BPJPSI_im", -100., 100.);
        addSU3BreakingParameter("delta_E2t_ccdd_BJPSIP_re", "E2t_ccsd_BJPSIP_re");
        addSU3BreakingParameter("delta_E2t_ccdd_BJPSIP_im", "E2t_ccsd_BJPSIP_im");
        addSU3BreakingParameter("delta_G2t_dcd_BJPSIP_re", "G2t_scd_BJPSIP_re");
        addSU3BreakingParameter("delta_G2t_dcd_BJPSIP_im", "G2t_scd_BJPSIP_im");
    }
    else if (channel == "Bdjpset8")
    {
        std::vector<std::string> params = {
            "delta_E2t_ccdd_BJPSIP_re",
            "delta_E2t_ccdd_BJPSIP_im",
            "dP4EW_ucd_BPJPSI_re",
            "dP4EW_ucd_BPJPSI_im",
            "EA2t_ccdd_BJPSIP_re",
            "EA2t_ccdd_BJPSIP_im",
            "delta_EA2t_ccsd_BJPSIP_re",
            "delta_EA2t_ccsd_BJPSIP_im",
            "EA2_ddcd_BPJPSI_re",
            "EA2_ddcd_BPJPSI_im",
            "delta_G2t_dcd_BJPSIP_re",
            "delta_G2t_dcd_BJPSIP_im",
            "G4t_cdd_BJPSIP_re",
            "G4t_cdd_BJPSIP_im",
            "delta_G4t_csd_BJPSIP_re",
            "delta_G4t_csd_BJPSIP_im"};
        channelParameters[channel] = params;

        addSU3BreakingParameter("delta_E2t_ccdd_BJPSIP_re", "E2t_ccsd_BJPSIP_re");
        addSU3BreakingParameter("delta_E2t_ccdd_BJPSIP_im", "E2t_ccsd_BJPSIP_im");
        addAmplitudeParameter("dP4EW_ucd_BPJPSI_re", -ewp_limit, ewp_limit);
        addAmplitudeParameter("dP4EW_ucd_BPJPSI_im", -ewp_limit, ewp_limit);
        addAmplitudeParameter("EA2t_ccdd_BJPSIP_re", -100., 100.);
        addAmplitudeParameter("EA2t_ccdd_BJPSIP_im", -100., 100.);
        addSU3BreakingParameter("delta_EA2t_ccsd_BJPSIP_re", "EA2t_ccdd_BJPSIP_re");
        addSU3BreakingParameter("delta_EA2t_ccsd_BJPSIP_im", "EA2t_ccdd_BJPSIP_im");
        addAmplitudeParameter("EA2_ddcd_BPJPSI_re", -100., 100.);
        addAmplitudeParameter("EA2_ddcd_BPJPSI_im", -100., 100.);
        addSU3BreakingParameter("delta_G2t_dcd_BJPSIP_re", "G2t_scd_BJPSIP_re");
        addSU3BreakingParameter("delta_G2t_dcd_BJPSIP_im", "G2t_scd_BJPSIP_im");
        addAmplitudeParameter("G4t_cdd_BJPSIP_re", -50., 50.);
        addAmplitudeParameter("G4t_cdd_BJPSIP_im", -50., 50.);
        addSU3BreakingParameter("delta_G4t_csd_BJPSIP_re", "G4t_cdd_BJPSIP_re");
        addSU3BreakingParameter("delta_G4t_csd_BJPSIP_im", "G4t_cdd_BJPSIP_im");
    }
    else if (channel == "Bdjpset1")
    {
        std::vector<std::string> params = {
            "delta_E2t_ccdd_BJPSIP_re",
            "delta_E2t_ccdd_BJPSIP_im",
            "dP4EW_ucd_BPJPSI_re",
            "dP4EW_ucd_BPJPSI_im",
            "EA2t_ccdd_BJPSIP_re",
            "EA2t_ccdd_BJPSIP_im",
            "delta_EA2t_ccsd_BJPSIP_re",
            "delta_EA2t_ccsd_BJPSIP_im",
            "EA2_ddcd_BPJPSI_re",
            "EA2_ddcd_BPJPSI_im",
            "delta_G2t_dcd_BJPSIP_re",
            "delta_G2t_dcd_BJPSIP_im",
            "G4t_cdd_BJPSIP_re",
            "G4t_cdd_BJPSIP_im",
            "delta_G4t_csd_BJPSIP_re",
            "delta_G4t_csd_BJPSIP_im"};
        channelParameters[channel] = params;

        addSU3BreakingParameter("delta_E2t_ccdd_BJPSIP_re", "E2t_ccsd_BJPSIP_re");
        addSU3BreakingParameter("delta_E2t_ccdd_BJPSIP_im", "E2t_ccsd_BJPSIP_im");
        addAmplitudeParameter("dP4EW_ucd_BPJPSI_re", -ewp_limit, ewp_limit);
        addAmplitudeParameter("dP4EW_ucd_BPJPSI_im", -ewp_limit, ewp_limit);
        addAmplitudeParameter("EA2t_ccdd_BJPSIP_re", -100., 100.);
        addAmplitudeParameter("EA2t_ccdd_BJPSIP_im", -100., 100.);
        addSU3BreakingParameter("delta_EA2t_ccsd_BJPSIP_re", "EA2t_ccdd_BJPSIP_re");
        addSU3BreakingParameter("delta_EA2t_ccsd_BJPSIP_im", "EA2t_ccdd_BJPSIP_im");
        addAmplitudeParameter("EA2_ddcd_BPJPSI_re", -100., 100.);
        addAmplitudeParameter("EA2_ddcd_BPJPSI_im", -100., 100.);
        addSU3BreakingParameter("delta_G2t_dcd_BJPSIP_re", "G2t_scd_BJPSIP_re");
        addSU3BreakingParameter("delta_G2t_dcd_BJPSIP_im", "G2t_scd_BJPSIP_im");
        addAmplitudeParameter("G4t_cdd_BJPSIP_re", -50., 50.);
        addAmplitudeParameter("G4t_cdd_BJPSIP_im", -50., 50.);
        addSU3BreakingParameter("delta_G4t_csd_BJPSIP_re", "G4t_cdd_BJPSIP_re");
        addSU3BreakingParameter("delta_G4t_csd_BJPSIP_im", "G4t_cdd_BJPSIP_im");
    }
    else if (channel = "Bpjpsikp")
    {
        std::vector<std::string> params = {
            "E2t_ccsd_BJPSIP_re",
            "E2t_ccsd_BJPSIP_im",
            "G2t_scd_BJPSIP_re",
            "G2t_scd_BJPSIP_im",
            "dp2EW_scu_BPJPSI_re",
            "dp2EW_scu_BPJPSI_im",
            "EA1_sdcd_BPJPSI_re",
            "EA1_sdcd_BPJPSI_im"};
        channelParameters[channel] = params;
        addAmplitudeParameter("E2t_ccsd_BJPSIP_re", 0., 100.);
        addAmplitudeParameter("E2t_ccsd_BJPSIP_im", 0., 0.);
        addAmplitudeParameter("G2t_scd_BJPSIP_re", -100., 100.);
        addAmplitudeParameter("G2t_scd_BJPSIP_im", -100., 100.);
        addAmplitudeParameter("dp2EW_scu_BPJPSI_re", -ewp_limit, ewp_limit);
        addAmplitudeParameter("dp2EW_scu_BPJPSI_im", -ewp_limit, ewp_limit);
        addAmplitudeParameter("EA1_sdcd_BPJPSI_re", -100., 100.);
        addAmplitudeParameter("EA1_sdcd_BPJPSI_im", -100., 100.);
    }
    else if (channel == "Bpjpsipp")
    {
        std::vector<std::string> params = {
            "delta_E2t_ccdd_BJPSIP_re",
            "delta_E2t_ccdd_BJPSIP_im",
            "delta_G2t_dcd_BJPSIP_re",
            "delta_G2t_dcd_BJPSIP_im",
            "delta_dp2EW_dcu_BPJPSI_re",
            "delta_dp2EW_dcu_BPJPSI_im",
            "delta_EA1_sdcd_BPJPSI_re",
            "delta_EA1_ddcd_BPJPSI_im"};
        channelParameters[channel] = params;
        addSU3BreakingParameter("delta_E2t_ccdd_BJPSIP_re", "E2t_ccsd_BJPSIP_re");
        addSU3BreakingParameter("delta_E2t_ccdd_BJPSIP_im", "E2t_ccsd_BJPSIP_im");
        addSU3BreakingParameter("delta_G2t_dcd_BJPSIP_re", "G2t_scd_BJPSIP_re");
        addSU3BreakingParameter("delta_G2t_dcd_BJPSIP_im", "G2t_scd_BJPSIP_im");
        addSU3BreakingParameter("delta_dp2EW_dcu_BPJPSI_re", "dp2EW_scu_BPJPSI_re");
        addSU3BreakingParameter("delta_dp2EW_dcu_BPJPSI_im", "dp2EW_scu_BPJPSI_im");
        addSU3BreakingParameter("delta_EA1_sdcd_BPJPSI_re", "EA1_sdcd_BPJPSI_re");
        addSU3BreakingParameter("delta_EA1_ddcd_BPJPSI_im", "EA1_sdcd_BPJPSI_im");
    }
    else if (channel == "Bsjpsip0")
    {
        std::vector<std::string> params = {
            "delta_dP4EW_ucs_BPJPSI_re",
            "delta_dP4EW_ucs_BPJPSI_im",
            "delta_EA2_ddcs_BPJPSI_re",
            "delta_EA2_ddcs_BPJPSI_im"};
        channelParameters[channel] = params;
        addSU3BreakingParameter("delta_dP4EW_ucs_BPJPSI_re", "dP4EW_ucd_BPJPSI_re");
        addSU3BreakingParameter("delta_dP4EW_ucs_BPJPSI_im", "dP4EW_ucd_BPJPSI_im");
        addSU3BreakingParameter("delta_EA2_ddcs_BPJPSI_re", "EA2_ddcd_BPJPSI_re");
        addSU3BreakingParameter("delta_EA2_ddcs_BPJPSI_im", "EA2_ddcd_BPJPSI_im");
    }
    else if (channel == "Bsjpsik0b")
    {
        std::vector<std::string> params = {
            "delta_E2t_ccdd_BJPSIP_re",
            "delta_E2t_ccdd_BJPSIP_im",
            "delta_E2t_ccds_BJPSIP_re",
            "delta_E2t_ccds_BJPSIP_im",
            "delta_G2t_dcd_BJPSIP_re",
            "delta_G2t_dcd_BJPSIP_im",
            "delta_G2t_dcs_BJPSIP_re",
            "delta_G2t_dcs_BJPSIP_im"};
        channelParameters[channel] = params;
        addSU3BreakingParameter("delta_E2t_ccdd_BJPSIP_re", "E2t_ccsd_BJPSIP_re");
        addSU3BreakingParameter("delta_E2t_ccdd_BJPSIP_im", "E2t_ccsd_BJPSIP_im");
        addSU3BreakingParameter("delta_E2t_ccds_BJPSIP_re", "delta_E2t_ccdd_BJPSIP_re");
        addSU3BreakingParameter("delta_E2t_ccds_BJPSIP_im", "delta_E2t_ccdd_BJPSIP_im");
        addSU3BreakingParameter("delta_G2t_dcd_BJPSIP_re", "G2t_scd_BJPSIP_re");
        addSU3BreakingParameter("delta_G2t_dcd_BJPSIP_im", "G2t_scd_BJPSIP_im");
        addSU3BreakingParameter("delta_G2t_dcs_BJPSIP_re", "delta_G2t_dcd_BJPSIP_re");
        addSU3BreakingParameter("delta_G2t_dcs_BJPSIP_im", "delta_G2t_dcd_BJPSIP_im");
    }
    else if (channel == "Bsjpsiet8")
    {
        std::vector<std::string> params = {
            "delta_E2t_ccss_BJPSIP_re",
            "delta_E2t_ccss_BJPSIP_im",
            "delta_P4EW_ucs_BPJPSI_re",
            "delta_P4EW_ucs_BPJPSI_im",
            "delta_EA2_ddcs_BPJPSI_re",
            "delta_EA2_ddcs_BPJPSI_im",
            "delta_G2t_scs_BJPSIP_re",
            "delta_G2t_scs_BJPSIP_im",
            "delta_G4t_cds_BJPSIP_re",
            "delta_G4t_cds_BJPSIP_im",
            "delta_G4t_css_BJPSIP_re",
            "delta_G4t_css_BJPSIP_im",
            "delta_EA2t_ccds_BJPSIP_re",
            "delta_EA2t_ccds_BJPSIP_im",
            "delta_EA2t_ccss_BJPSIP_re",
            "delta_EA2t_ccss_BJPSIP_im"};
        channelParameters[channel] = params;
        addSU3BreakingParameter("delta_E2t_ccss_BJPSIP_re", "E2t_ccsd_BJPSIP_re");
        addSU3BreakingParameter("delta_E2t_ccss_BJPSIP_im", "E2t_ccsd_BJPSIP_im");
        addSU3BreakingParameter("delta_P4EW_ucs_BPJPSI_re", "dP4EW_ucd_BPJPSI_re");
        addSU3BreakingParameter("delta_P4EW_ucs_BPJPSI_im", "dP4EW_ucd_BPJPSI_im");
        addSU3BreakingParameter("delta_EA2_ddcs_BPJPSI_re", "EA2_ddcd_BPJPSI_re");
        addSU3BreakingParameter("delta_EA2_ddcs_BPJPSI_im", "EA2_ddcd_BPJPSI_im");
        addSU3BreakingParameter("delta_G2t_scs_BJPSIP_re", "G2t_scd_BJPSIP_re");
        addSU3BreakingParameter("delta_G2t_scs_BJPSIP_im", "G2t_scd_BJPSIP_im");
        addSU3BreakingParameter("delta_G4t_cds_BJPSIP_re", "G4t_cdd_BJPSIP_re");
        addSU3BreakingParameter("delta_G4t_cds_BJPSIP_im", "G4t_cdd_BJPSIP_im");
        addSU3BreakingParameter("delta_G4t_css_BJPSIP_re", "delta_G4t_cds_BJPSIP_re");
        addSU3BreakingParameter("delta_G4t_css_BJPSIP_im", "delta_G4t_cds_BJPSIP_im");
        addSU3BreakingParameter("delta_EA2t_ccds_BJPSIP_re", "EA2t_ccdd_BJPSIP_re");
        addSU3BreakingParameter("delta_EA2t_ccds_BJPSIP_im", "EA2t_ccdd_BJPSIP_im");
        addSU3BreakingParameter("delta_EA2t_ccss_BJPSIP_re", "delta_EA2t_ccds_BJPSIP_re");
        addSU3BreakingParameter("delta_EA2t_ccss_BJPSIP_im", "delta_EA2t_ccds_BJPSIP_im");
    }
    else if (channel == "Bsjpsiet1")
    {
        std::vector<std::string> params = {
            "delta_E2t_ccss_BJPSIP_re",
            "delta_E2t_ccss_BJPSIP_im",
            "delta_P4EW_ucs_BPJPSI_re",
            "delta_P4EW_ucs_BPJPSI_im",
            "delta_EA2_ddcs_BPJPSI_re",
            "delta_EA2_ddcs_BPJPSI_im",
            "delta_G2t_scs_BJPSIP_re",
            "delta_G2t_scs_BJPSIP_im",
            "delta_G4t_cds_BJPSIP_re",
            "delta_G4t_cds_BJPSIP_im",
            "delta_G4t_css_BJPSIP_re",
            "delta_G4t_css_BJPSIP_im",
            "delta_EA2t_ccds_BJPSIP_re",
            "delta_EA2t_ccds_BJPSIP_im",
            "delta_EA2t_ccss_BJPSIP_re",
            "delta_EA2t_ccss_BJPSIP_im"};
        channelParameters[channel] = params;
        addSU3BreakingParameter("delta_E2t_ccss_BJPSIP_re", "E2t_ccsd_BJPSIP_re");
        addSU3BreakingParameter("delta_E2t_ccss_BJPSIP_im", "E2t_ccsd_BJPSIP_im");
        addSU3BreakingParameter("delta_P4EW_ucs_BPJPSI_re", "dP4EW_ucd_BPJPSI_re");
        addSU3BreakingParameter("delta_P4EW_ucs_BPJPSI_im", "dP4EW_ucd_BPJPSI_im");
        addSU3BreakingParameter("delta_EA2_ddcs_BPJPSI_re", "EA2_ddcd_BPJPSI_re");
        addSU3BreakingParameter("delta_EA2_ddcs_BPJPSI_im", "EA2_ddcd_BPJPSI_im");
        addSU3BreakingParameter("delta_G2t_scs_BJPSIP_re", "G2t_scd_BJPSIP_re");
        addSU3BreakingParameter("delta_G2t_scs_BJPSIP_im", "G2t_scd_BJPSIP_im");
        addSU3BreakingParameter("delta_G4t_cds_BJPSIP_re", "G4t_cdd_BJPSIP_re");
        addSU3BreakingParameter("delta_G4t_cds_BJPSIP_im", "G4t_cdd_BJPSIP_im");
        addSU3BreakingParameter("delta_G4t_css_BJPSIP_re", "delta_G4t_cds_BJPSIP_re");
        addSU3BreakingParameter("delta_G4t_css_BJPSIP_im", "delta_G4t_cds_BJPSIP_im");
        addSU3BreakingParameter("delta_EA2t_ccds_BJPSIP_re", "EA2t_ccdd_BJPSIP_re");
        addSU3BreakingParameter("delta_EA2t_ccds_BJPSIP_im", "EA2t_ccdd_BJPSIP_im");
        addSU3BreakingParameter("delta_EA2t_ccss_BJPSIP_re", "delta_EA2t_ccds_BJPSIP_re");
        addSU3BreakingParameter("delta_EA2t_ccss_BJPSIP_im", "delta_EA2t_ccds_BJPSIP_im");
    }
    // Vector channels with helicity amplitudes
    else if (channel == "Bsjpsiph")
    {
        std::vector<std::string> params = {
            "E2t_ccss_BJPSIV_0_re",
            "E2t_ccss_BJPSIV_0_im",
            "G2t_scs_BJPSIV_0_re",
            "G2t_scs_BJPSIV_0_im",
            "EA2t_ccss_BJPSIV_0_re",
            "EA2t_ccss_BJPSIV_0_im",
            "G4t_css_BJPSIV_0_re",
            "G4t_css_BJPSIV_0_im",
            "E2t_ccss_BJPSIV_paral_re",
            "E2t_ccss_BJPSIV_paral_im",
            "G2t_scs_BJPSIV_paral_re",
            "G2t_scs_BJPSIV_paral_im",
            "EA2t_ccss_BJPSIV_paral_re",
            "EA2t_ccss_BJPSIV_paral_im",
            "G4t_css_BJPSIV_paral_re",
            "G4t_css_BJPSIV_paral_im",
            "E2t_ccss_BJPSIV_perp_re",
            "E2t_ccss_BJPSIV_perp_im",
            "G2t_scs_BJPSIV_perp_re",
            "G2t_scs_BJPSIV_perp_im",
            "EA2t_ccss_BJPSIV_perp_re",
            "EA2t_ccss_BJPSIV_perp_im",
            "G4t_css_BJPSIV_perp_re",
            "G4t_css_BJPSIV_perp_im"};
        channelParameters[channel] = params;

        addAmplitudeParameter("E2t_ccss_BJPSIV_re", 0., 100., true);
        addAmplitudeParameter("E2t_ccss_BJPSIV_im", 0., 0., true);
        addAmplitudeParameter("G2t_scs_BJPSIV_re", -100., 100., true);
        addAmplitudeParameter("G2t_scs_BJPSIV_im", -100., 100., true);
        addAmplitudeParameter("EA2t_ccss_BJPSIV_re", -100., 100., true);
        addAmplitudeParameter("EA2t_ccss_BJPSIV_im", -100., 100., true);
        addAmplitudeParameter("G4t_css_BJPSIV_re", -50., 50., true);
        addAmplitudeParameter("G4t_css_BJPSIV_im", -50., 50., true);
    }
    else if (channel == "Bsjpsiom")
    {
        std::vector<std::string> params = {
            "delta_EA2t_ccds_BJPSIV_0_re",
            "delta_EA2t_ccds_BJPSIV_0_im",
            "delta_G4t_cds_BJPSIV_0_re",
            "delta_G4t_cds_BJPSIV_0_im",
            "dp4EW_ucs_BVJPSI_0_re",
            "dp4EW_ucs_BVJPSI_0_im",
            "EA2_ddcs_BVJPSI_0_re",
            "EA2_ddcs_BVJPSI_0_im",
            "delta_EA2t_ccds_BJPSIV_paral_re",
            "delta_EA2t_ccds_BJPSIV_paral_im",
            "delta_G4t_cds_BJPSIV_paral_re",
            "delta_G4t_cds_BJPSIV_paral_im",
            "dp4EW_ucs_BVJPSI_paral_re",
            "dp4EW_ucs_BVJPSI_paral_im",
            "EA2_ddcs_BVJPSI_paral_re",
            "EA2_ddcs_BVJPSI_paral_im",
            "delta_EA2t_ccds_BJPSIV_perp_re",
            "delta_EA2t_ccds_BJPSIV_perp_im",
            "delta_G4t_cds_BJPSIV_perp_re",
            "delta_G4t_cds_BJPSIV_perp_im",
            "dp4EW_ucs_BVJPSI_perp_re",
            "dp4EW_ucs_BVJPSI_perp_im",
            "EA2_ddcs_BVJPSI_perp_re",
            "EA2_ddcs_BVJPSI_perp_im"};
        channelParameters[channel] = params;
        addSU3BreakingParameter("delta_EA2t_ccds_BJPSIV_re", "EA2t_ccss_BJPSIV_re", true);
        addSU3BreakingParameter("delta_EA2t_ccds_BJPSIV_im", "EA2t_ccss_BJPSIV_im", true);
        addSU3BreakingParameter("delta_G4t_cds_BJPSIV_re", "G4t_css_BJPSIV_re", true);
        addSU3BreakingParameter("delta_G4t_cds_BJPSIV_im", "G4t_css_BJPSIV_im", true);
        addAmplitudeParameter("dp4EW_ucs_BVJPSI_re", -ewp_limit, ewp_limit, true);
        addAmplitudeParameter("dp4EW_ucs_BVJPSI_im", -ewp_limit, ewp_limit, true);
        addAmplitudeParameter("EA2_ddcs_BVJPSI_re", -100., 100., true);
        addAmplitudeParameter("EA2_ddcs_BVJPSI_im", -100., 100., true);
    }
    else if (channel == "Bsjpsikbst")
    {
        std::vector<std::string> params = {
            "delta_E2t_ccds_BJPSIV_0_re",
            "delta_E2t_ccds_BJPSIV_0_im",
            "delta_G2t_dcs_BJPSIV_0_re",
            "delta_G2t_dcs_BJPSIV_0_im",
            "delta_E2t_ccds_BJPSIV_paral_re",
            "delta_E2t_ccds_BJPSIV_paral_im",
            "delta_G2t_dcs_BJPSIV_paral_re",
            "delta_G2t_dcs_BJPSIV_paral_im",
            "delta_E2t_ccds_BJPSIV_perp_re",
            "delta_E2t_ccds_BJPSIV_perp_im",
            "delta_G2t_dcs_BJPSIV_perp_re",
            "delta_G2t_dcs_BJPSIV_perp_im"};
        channelParameters[channel] = params;

        addSU3BreakingParameter("delta_E2t_ccds_BJPSIV_re", "E2t_ccss_BJPSIV_re", true);
        addSU3BreakingParameter("delta_E2t_ccds_BJPSIV_im", "E2t_ccss_BJPSIV_im", true);
        addSU3BreakingParameter("delta_G2t_dcs_BJPSIV_re", "G2t_scs_BJPSIV_re", true);
        addSU3BreakingParameter("delta_G2t_dcs_BJPSIV_im", "G2t_scs_BJPSIV_im", true);
    }
    else if (channel == "Bsjpsirho0")
    {
        std::vector<std::string> params = {
            "dP4EW_ucs_BVJPSI_0_re",
            "dP4EW_ucs_BVJPSI_0_im",
            "EA2_ddcs_BVJPSI_0_re",
            "EA2_ddcs_BVJPSI_0_im",
            "dP4EW_ucs_BVJPSI_paral_re",
            "dP4EW_ucs_BVJPSI_paral_im",
            "EA2_ddcs_BVJPSI_paral_re",
            "EA2_ddcs_BVJPSI_paral_im",
            "dP4EW_ucs_BVJPSI_perp_re",
            "dP4EW_ucs_BVJPSI_perp_im",
            "EA2_ddcs_BVJPSI_perp_re",
            "EA2_ddcs_BVJPSI_perp_im"};
        channelParameters[channel] = params;

        addAmplitudeParameter("dp4EW_ucs_BVJPSI_re", -ewp_limit, ewp_limit, true);
        addAmplitudeParameter("dp4EW_ucs_BVJPSI_im", -ewp_limit, ewp_limit, true);
        addAmplitudeParameter("EA2_ddcs_BVJPSI_re", -100., 100., true);
        addAmplitudeParameter("EA2_ddcs_BVJPSI_im", -100., 100., true);
    }
    else if (channel == "Bdjpsiom")
    {
        std::vector<std::string> params = {
            "delta_E2t_ccds_BJPSIV_0_re",
            "delta_E2t_ccds_BJPSIV_0_im",
            "delta_E2t_ccdd_BJPSIV_0_re",
            "delta_E2t_ccdd_BJPSIV_0_im",
            "delta_G2t_dcs_BJPSIV_0_re",
            "delta_G2t_dcs_BJPSIV_0_im",
            "delta_G2t_dcd_BJPSIV_0_re",
            "delta_G2t_dcd_BJPSIV_0_im",
            "delta_dp4EW_ucd_BVJPSI_0_re",
            "delta_dp4EW_ucd_BVJPSI_0_im",
            "delta_EA2t_ccds_BJPSIV_0_re",
            "delta_EA2t_ccds_BJPSIV_0_im",
            "delta_EA2t_ccdd_BJPSIV_0_re",
            "delta_EA2t_ccdd_BJPSIV_0_im",
            "delta_EA2_ddcd_BVJPSI_0_re",
            "delta_EA2_ddcd_BVJPSI_0_im",
            "delta_G4t_cds_BJPSIV_0_re",
            "delta_G4t_cds_BJPSIV_0_im",
            "delta_G4t_cdd_BJPSIV_0_re",
            "delta_G4t_cdd_BJPSIV_0_im",
            "delta_E2t_ccds_BJPSIV_paral_re",
            "delta_E2t_ccds_BJPSIV_paral_im",
            "delta_E2t_ccdd_BJPSIV_paral_re",
            "delta_E2t_ccdd_BJPSIV_paral_im",
            "delta_G2t_dcs_BJPSIV_paral_re",
            "delta_G2t_dcs_BJPSIV_paral_im",
            "delta_G2t_dcd_BJPSIV_paral_re",
            "delta_G2t_dcd_BJPSIV_paral_im",
            "delta_dp4EW_ucd_BVJPSI_paral_re",
            "delta_dp4EW_ucd_BVJPSI_paral_im",
            "delta_EA2t_ccds_BJPSIV_paral_re",
            "delta_EA2t_ccds_BJPSIV_paral_im",
            "delta_EA2t_ccdd_BJPSIV_paral_re",
            "delta_EA2t_ccdd_BJPSIV_paral_im",
            "delta_EA2_ddcd_BVJPSI_paral_re",
            "delta_EA2_ddcd_BVJPSI_paral_im",
            "delta_G4t_cds_BJPSIV_paral_re",
            "delta_G4t_cds_BJPSIV_paral_im",
            "delta_G4t_cdd_BJPSIV_paral_re",
            "delta_G4t_cdd_BJPSIV_paral_im",
            "delta_E2t_ccds_BJPSIV_perp_re",
            "delta_E2t_ccds_BJPSIV_perp_im",
            "delta_E2t_ccdd_BJPSIV_perp_re",
            "delta_E2t_ccdd_BJPSIV_perp_im",
            "delta_G2t_dcs_BJPSIV_perp_re",
            "delta_G2t_dcs_BJPSIV_perp_im",
            "delta_G2t_dcd_BJPSIV_perp_re",
            "delta_G2t_dcd_BJPSIV_perp_im",
            "delta_dp4EW_ucd_BVJPSI_perp_re",
            "delta_dp4EW_ucd_BVJPSI_perp_im",
            "delta_EA2t_ccds_BJPSIV_perp_re",
            "delta_EA2t_ccds_BJPSIV_perp_im",
            "delta_EA2t_ccdd_BJPSIV_perp_re",
            "delta_EA2t_ccdd_BJPSIV_perp_im",
            "delta_EA2_ddcd_BVJPSI_perp_re",
            "delta_EA2_ddcd_BVJPSI_perp_im",
            "delta_G4t_cds_BJPSIV_perp_re",
            "delta_G4t_cds_BJPSIV_perp_im",
            "delta_G4t_cdd_BJPSIV_perp_re",
            "delta_G4t_cdd_BJPSIV_perp_im"};

        channelParameters[channel] = params;

        addSU3BreakingParameter("delta_E2t_ccds_BJPSIV_re", "E2t_ccss_BJPSIV_re", true);
        addSU3BreakingParameter("delta_E2t_ccds_BJPSIV_im", "E2t_ccss_BJPSIV_im", true);
        addSU3BreakingParameter("delta_E2t_ccdd_BJPSIV_re", "E2t_ccds_BJPSIV_re", true);
        addSU3BreakingParameter("delta_E2t_ccdd_BJPSIV_im", "E2t_ccds_BJPSIV_im", true);
        addSU3BreakingParameter("delta_G2t_dcs_BJPSIV_re", "G2t_scs_BJPSIV_re", true);
        addSU3BreakingParameter("delta_G2t_dcs_BJPSIV_im", "G2t_scs_BJPSIV_im", true);
        addSU3BreakingParameter("delta_G2t_dcd_BJPSIV_re", "G2t_dcs_BJPSIV_re", true);
        addSU3BreakingParameter("delta_G2t_dcd_BJPSIV_im", "G2t_dcs_BJPSIV_im", true);
        addSU3BreakingParameter("delta_dp4EW_ucd_BVJPSI_re", "dp4EW_ucs_BVJPSI_re", true);
        addSU3BreakingParameter("delta_dp4EW_ucd_BVJPSI_im", "dp4EW_ucs_BVJPSI_im", true);
        addSU3BreakingParameter("delta_EA2t_ccds_BJPSIV_re", "EA2t_ccss_BJPSIV_re", true);
        addSU3BreakingParameter("delta_EA2t_ccds_BJPSIV_im", "EA2t_ccss_BJPSIV_im", true);
        addSU3BreakingParameter("delta_EA2t_ccdd_BJPSIV_re", "EA2t_ccds_BJPSIV_re", true);
        addSU3BreakingParameter("delta_EA2t_ccdd_BJPSIV_im", "EA2t_ccds_BJPSIV_im", true);
        addSU3BreakingParameter("delta_EA2_ddcd_BVJPSI_re", "EA2_ddcs_BVJPSI_re", true);
        addSU3BreakingParameter("delta_EA2_ddcd_BVJPSI_im", "EA2_ddcs_BVJPSI_im", true);
        addSU3BreakingParameter("delta_G4t_cds_BJPSIV_re", "G4t_css_BJPSIV_re", true);
        addSU3BreakingParameter("delta_G4t_cds_BJPSIV_im", "G4t_css_BJPSIV_im", true);
        addSU3BreakingParameter("delta_G4t_cdd_BJPSIV_re", "G4t_cds_BJPSIV_re", true);
        addSU3BreakingParameter("delta_G4t_cdd_BJPSIV_im", "G4t_cds_BJPSIV_im", true);
    }
    else if (channel == "Bdjpsikst")
    {
        std::vector<std::string> params = {
            "delta_E2t_ccsd_BJPSIV_0_re",
            "delta_E2t_ccsd_BJPSIV_0_im",
            "delta_G2t_scd_BJPSIV_0_re",
            "delta_G2t_scd_BJPSIV_0_im",
            "delta_E2t_ccsd_BJPSIV_paral_re",
            "delta_E2t_ccsd_BJPSIV_paral_im",
            "delta_G2t_scd_BJPSIV_paral_re",
            "delta_G2t_scd_BJPSIV_paral_im",
            "delta_E2t_ccsd_BJPSIV_perp_re",
            "delta_E2t_ccsd_BJPSIV_perp_im",
            "delta_G2t_scd_BJPSIV_perp_re",
            "delta_G2t_scd_BJPSIV_perp_im"};
        channelParameters[channel] = params;
        addSU3BreakingParameter("delta_E2t_ccsd_BJPSIV_re", "E2t_ccss_BJPSIV_re", true);
        addSU3BreakingParameter("delta_E2t_ccsd_BJPSIV_im", "E2t_ccss_BJPSIV_im", true);
        addSU3BreakingParameter("delta_G2t_scd_BJPSIV_re", "G2t_scs_BJPSIV_re", true);
        addSU3BreakingParameter("delta_G2t_scd_BJPSIV_im", "G2t_scs_BJPSIV_im", true);
    }
    else if (channel == "Bdjpsirho0")
    {
        std::vector<std::string> params = {
            "delta_E2t_ccsd_BJPSIV_0_re",
            "delta_E2t_ccsd_BJPSIV_0_im",
            "delta_E2t_ccdd_BJPSIV_0_re",
            "delta_E2t_ccdd_BJPSIV_0_im",
            "delta_P4EW_ucd_BVJPSI_0_re",
            "delta_P4EW_ucd_BVJPSI_0_im",
            "delta_EA2_ddcd_BVJPSI_0_re",
            "delta_EA2_ddcd_BVJPSI_0_im",
            "delta_G2t_dcs_BJPSIV_0_re",
            "delta_G2t_dcs_BJPSIV_0_im",
            "delta_G2t_dcd_BJPSIV_0_re",
            "delta_G2t_dcd_BJPSIV_0_im",
            "delta_E2t_ccsd_BJPSIV_paral_re",
            "delta_E2t_ccsd_BJPSIV_paral_im",
            "delta_E2t_ccdd_BJPSIV_paral_re",
            "delta_E2t_ccdd_BJPSIV_paral_im",
            "delta_P4EW_ucd_BVJPSI_paral_re",
            "delta_P4EW_ucd_BVJPSI_paral_im",
            "delta_EA2_ddcd_BVJPSI_paral_re",
            "delta_EA2_ddcd_BVJPSI_paral_im",
            "delta_G2t_dcs_BJPSIV_paral_re",
            "delta_G2t_dcs_BJPSIV_paral_im",
            "delta_G2t_dcd_BJPSIV_paral_re",
            "delta_G2t_dcd_BJPSIV_paral_im",
            "delta_E2t_ccsd_BJPSIV_perp_re",
            "delta_E2t_ccsd_BJPSIV_perp_im",
            "delta_E2t_ccdd_BJPSIV_perp_re",
            "delta_E2t_ccdd_BJPSIV_perp_im",
            "delta_P4EW_ucd_BVJPSI_perp_re",
            "delta_P4EW_ucd_BVJPSI_perp_im",
            "delta_EA2_ddcd_BVJPSI_perp_re",
            "delta_EA2_ddcd_BVJPSI_perp_im",
            "delta_G2t_dcs_BJPSIV_perp_re",
            "delta_G2t_dcs_BJPSIV_perp_im",
            "delta_G2t_dcd_BJPSIV_perp_re",
            "delta_G2t_dcd_BJPSIV_perp_im"};
        channelParameters[channel] = params;

        addSU3BreakingParameter("delta_E2t_ccsd_BJPSIV_re", "E2t_ccss_BJPSIV_re", true);
        addSU3BreakingParameter("delta_E2t_ccsd_BJPSIV_im", "E2t_ccss_BJPSIV_im", true);
        addSU3BreakingParameter("delta_E2t_ccdd_BJPSIV_re", "E2t_ccsd_BJPSIV_re", true);
        addSU3BreakingParameter("delta_E2t_ccdd_BJPSIV_im", "E2t_ccsd_BJPSIV_im", true);
        addSU3BreakingParameter("delta_P4EW_ucd_BVJPSI_re", "dP4EW_ucs_BVJPSI_re", true);
        addSU3BreakingParameter("delta_P4EW_ucd_BVJPSI_im", "dP4EW_ucs_BVJPSI_im", true);
        addSU3BreakingParameter("delta_EA2_ddcd_BVJPSI_re", "EA2_ddcs_BVJPSI_re", true);
        addSU3BreakingParameter("delta_EA2_ddcd_BVJPSI_im", "EA2_ddcs_BVJPSI_im", true);
        addSU3BreakingParameter("delta_G2t_dcs_BJPSIV_re", "G2t_scs_BJPSIV_re", true);
        addSU3BreakingParameter("delta_G2t_dcs_BJPSIV_im", "G2t_scs_BJPSIV_im", true);
        addSU3BreakingParameter("delta_G2t_dcd_BJPSIV_re", "G2t_dcs_BJPSIV_re", true);
        addSU3BreakingParameter("delta_G2t_dcd_BJPSIV_im", "G2t_dcs_BJPSIV_im", true);
    }
    else if (channel == "Bdjpsiph")
    {
        std::vector<std::string> params = {
            "delta_EA2t_ccsd_BJPSIV_0_re",
            "delta_EA2t_ccsd_BJPSIV_0_im",
            "delta_G4t_csd_BJPSIV_0_re",
            "delta_G4t_csd_BJPSIV_0_im",
            "delta_EA2t_ccsd_BJPSIV_paral_re",
            "delta_EA2t_ccsd_BJPSIV_paral_im",
            "delta_G4t_csd_BJPSIV_paral_re",
            "delta_G4t_csd_BJPSIV_paral_im",
            "delta_EA2t_ccsd_BJPSIV_perp_re",
            "delta_EA2t_ccsd_BJPSIV_perp_im",
            "delta_G4t_csd_BJPSIV_perp_re",
            "delta_G4t_csd_BJPSIV_perp_im"};
        channelParameters[channel] = params;

        addSU3BreakingParameter("delta_EA2t_ccsd_BJPSIV_re", "EA2t_ccss_BJPSIV_re", true);
        addSU3BreakingParameter("delta_EA2t_ccsd_BJPSIV_im", "EA2t_ccss_BJPSIV_im", true);
        addSU3BreakingParameter("delta_G4t_csd_BJPSIV_re", "G4t_css_BJPSIV_re", true);
        addSU3BreakingParameter("delta_G4t_csd_BJPSIV_im", "G4t_css_BJPSIV_im", true);
    }
    else if (channel == "Bpjpsikst")
    {
        std::vector<std::string> params = {
            "delta_E2t_ccsd_BJPSIV_0_re",
            "delta_E2t_ccsd_BJPSIV_0_im",
            "dP2EW_scu_BJPSIV_0_re",
            "dP2EW_scu_BJPSIV_0_im",
            "EA1_sdcd_BVJPSI_0_re",
            "EA1_sdcd_BVJPSI_0_im",
            "delta_G2t_scd_BJPSIV_0_re",
            "delta_G2t_scd_BJPSIV_0_im",
            "delta_E2t_ccsd_BJPSIV_paral_re",
            "delta_E2t_ccsd_BJPSIV_paral_im",
            "dP2EW_scu_BJPSIV_paral_re",
            "dP2EW_scu_BJPSIV_paral_im",
            "EA1_sdcd_BVJPSI_paral_re",
            "EA1_sdcd_BVJPSI_paral_im",
            "delta_G2t_scd_BJPSIV_paral_re",
            "delta_G2t_scd_BJPSIV_paral_im",
            "delta_E2t_ccsd_BJPSIV_perp_re",
            "delta_E2t_ccsd_BJPSIV_perp_im",
            "dP2EW_scu_BJPSIV_perp_re",
            "dP2EW_scu_BJPSIV_perp_im",
            "EA1_sdcd_BVJPSI_perp_re",
            "EA1_sdcd_BVJPSI_perp_im",
            "delta_G2t_scd_BJPSIV_perp_re",
            "delta_G2t_scd_BJPSIV_perp_im"};

        channelParameters[channel] = params;

        addSU3BreakingParameter("delta_E2t_ccsd_BJPSIV_re", "E2t_ccss_BJPSIV_re", true);
        addSU3BreakingParameter("delta_E2t_ccsd_BJPSIV_im", "E2t_ccss_BJPSIV_im", true);
        addAmplitudeParameter("dP2EW_scu_BJPSIV_re", -ewp_limit, ewp_limit, true);
        addAmplitudeParameter("dP2EW_scu_BJPSIV_im", -ewp_limit, ewp_limit, true);
        addAmplitudeParameter("EA1_sdcd_BVJPSI_re", -100., 100., true);
        addAmplitudeParameter("EA1_sdcd_BVJPSI_im", -100., 100., true);
        addSU3BreakingParameter("delta_G2t_scd_BJPSIV_re", "G2t_scs_BJPSIV_re", true);
        addSU3BreakingParameter("delta_G2t_scd_BJPSIV_im", "G2t_scs_BJPSIV_im", true);
    }
    else if (channel == "Bpjpsirhop")
    {
        std::vector<std::string> params = {
            "delta_E2t_ccsd_BJPSIV_0_re",
            "delta_E2t_ccsd_BJPSIV_0_im",
            "delta_E2t_ccdd_BJPSIV_0_re",
            "delta_E2t_ccdd_BJPSIV_0_im",
            "delta_dP2EW_dcu_BJPSIV_0_re",
            "delta_dP2EW_dcu_BJPSIV_0_im",
            "delta_EA1_ddcd_BVJPSI_0_re",
            "delta_EA1_ddcd_BVJPSI_0_im",
            "delta_G2t_scd_BJPSIV_0_re",
            "delta_G2t_scd_BJPSIV_0_im",
            "delta_G2t_dcd_BJPSIV_0_re",
            "delta_G2t_dcd_BJPSIV_0_im",
            "delta_E2t_ccsd_BJPSIV_paral_re",
            "delta_E2t_ccsd_BJPSIV_paral_im",
            "delta_E2t_ccdd_BJPSIV_paral_re",
            "delta_E2t_ccdd_BJPSIV_paral_im",
            "delta_dP2EW_dcu_BJPSIV_paral_re",
            "delta_dP2EW_dcu_BJPSIV_paral_im",
            "delta_EA1_ddcd_BVJPSI_paral_re",
            "delta_EA1_ddcd_BVJPSI_paral_im",
            "delta_G2t_scd_BJPSIV_paral_re",
            "delta_G2t_scd_BJPSIV_paral_im",
            "delta_G2t_dcd_BJPSIV_paral_re",
            "delta_G2t_dcd_BJPSIV_paral_im",
            "delta_E2t_ccsd_BJPSIV_perp_re",
            "delta_E2t_ccsd_BJPSIV_perp_im",
            "delta_E2t_ccdd_BJPSIV_perp_re",
            "delta_E2t_ccdd_BJPSIV_perp_im",
            "delta_dP2EW_dcu_BJPSIV_perp_re",
            "delta_dP2EW_dcu_BJPSIV_perp_im",
            "delta_EA1_ddcd_BVJPSI_perp_re",
            "delta_EA1_ddcd_BVJPSI_perp_im",
            "delta_G2t_scd_BJPSIV_perp_re",
            "delta_G2t_scd_BJPSIV_perp_im",
            "delta_G2t_dcd_BJPSIV_perp_re",
            "delta_G2t_dcd_BJPSIV_perp_im"};

        channelParameters[channel] = params;

        addSU3BreakingParameter("delta_E2t_ccsd_BJPSIV_re", "E2t_ccss_BJPSIV_re", true);
        addSU3BreakingParameter("delta_E2t_ccsd_BJPSIV_im", "E2t_ccss_BJPSIV_im", true);
        addSU3BreakingParameter("delta_E2t_ccdd_BJPSIV_re", "E2t_ccsd_BJPSIV_re", true);
        addSU3BreakingParameter("delta_E2t_ccdd_BJPSIV_im", "E2t_ccsd_BJPSIV_im", true);
        addSU3BreakingParameter("delta_dP2EW_dcu_BJPSIV_re", "dP2EW_scu_BJPSIV_re", true);
        addSU3BreakingParameter("delta_dP2EW_dcu_BJPSIV_im", "dP2EW_scu_BJPSIV_im", true);
        addSU3BreakingParameter("delta_EA1_ddcd_BVJPSI_re", "EA1_sdcd_BVJPSI_re", true);
        addSU3BreakingParameter("delta_EA1_ddcd_BVJPSI_im", "EA1_sdcd_BVJPSI_im", true);
        addSU3BreakingParameter("delta_G2t_scd_BJPSIV_re", "G2t_scs_BJPSIV_re", true);
        addSU3BreakingParameter("delta_G2t_scd_BJPSIV_im", "G2t_scs_BJPSIV_im", true);
        addSU3BreakingParameter("delta_G2t_dcd_BJPSIV_re", "G2t_dcs_BJPSIV_re", true);
        addSU3BreakingParameter("delta_G2t_dcd_BJPSIV_im", "G2t_dcs_BJPSIV_im", true);
    }
    else if (channel == "Bsdpsdms")
    {
        // Bs -> Ds+ Ds- parameters (base channel for s->c)
        // b → c(c̄s), spectator s
        std::vector<std::string> params = {
            "E1t_sccs_BDDb_re", "E1t_sccs_BDDb_im",
            "A2t_cscs_BDbD_re", "A2t_cscs_BDbD_im",
            "G1t_scs_BDDb_re", "G1t_scs_BDDb_im",
            "G3t_css_BDDb_re", "G3t_css_BDDb_im"};
        channelParameters[channel] = params;

        addAmplitudeParameter("E1t_sccs_BDDb_re", -100., 100.);
        addAmplitudeParameter("E1t_sccs_BDDb_im", 0., 0.);
        addAmplitudeParameter("A2t_cscs_BDbD_re", -100., 100.);
        addAmplitudeParameter("A2t_cscs_BDbD_im", -100., 100.);
        addAmplitudeParameter("G1t_scs_BDDb_re", -100., 100.);
        addAmplitudeParameter("G1t_scs_BDDb_im", -100., 100.);
        addAmplitudeParameter("G3t_css_BDDb_re", -100., 100.);
        addAmplitudeParameter("G3t_css_BDDb_im", -100., 100.);
    }
    else if (channel == "Bsdpdms")
    {
        // b → c(c̄d), spectator s
        std::vector<std::string> params = {
            "delta_E1t_dccs_BDDb_re", "delta_E1t_dccs_BDDb_im",
            "delta_G1t_dcs_BDDb_re", "delta_G1t_dcs_BDDb_im"};
        channelParameters[channel] = params;

        AddSU3BreakingParameter("delta_E1t_dccs_BDDb_re", "E1t_sccs_BDDb_re");
        AddSU3BreakingParameter("delta_E1t_dccs_BDDb_im", "E1t_sccs_BDDb_im");
        AddSU3BreakingParameter("delta_G1t_dcs_BDDb_re", "G1t_scs_BDDb_re");
        AddSU3BreakingParameter("delta_G1t_dcs_BDDb_im", "G1t_scs_BDDb_im");
    }
    else if (channel == "Bsdpdm")
    {
        // b → c(c̄s), spectator s
        std::vector<std::string> params = {
            "delta_A2t_cdcs_BDbD_re", "delta_A2t_cdcs_BDbD_im",
            "delta_G3_cds_BDDb_re", "delta_G3_cds_BDDb_im"};
        channelParameters[channel] = params;

        AddSU3BreakingParameter("delta_A2t_cdcs_BDbD_re", "A2t_cscs_BDbD_re");
        AddSU3BreakingParameter("delta_A2t_cdcs_BDbD_im", "A2t_cscs_BDbD_im");
        AddSU3BreakingParameter("delta_G3_cds_BDDb_re", "G3t_css_BDDb_re");
        AddSU3BreakingParameter("delta_G3_cds_BDDb_im", "G3t_css_BDDb_im");
    }
    else if (channel == "Bsd0d0b")
    {
        // b → c(c̄s), spectator s
        std::vector<std::string> params = {
            "delta_A2t_cdcs_BDbD_re", "delta_A2t_cdcs_BDbD_im",
            "dP3EW_ucs_BDbD_re", "dP3EW_ucs_BDbD_im",
            "A2_dcds_BDDb_re", "A2_dcds_BDDb_im",
            "delta_G3t_cds_BDDb_re", "delta_G3t_cds_BDDb_im"};
        channelParameters[channel] = params;

        AddSU3BreakingParameter("delta_A2t_cdcs_BDbD_re", "A2t_cscs_BDbD_re");
        AddSU3BreakingParameter("delta_A2t_cdcs_BDbD_im", "A2t_cscs_BDbD_im");
        addAmplitudeParameter("dP3EW_ucs_BDbD_re", -ewp_limit, ewp_limit);
        addAmplitudeParameter("dP3EW_ucs_BDbD_im", -ewp_limit, ewp_limit);
        addAmplitudeParameter("A2_dcds_BDDb_re", -100., 100.);
        addAmplitudeParameter("A2_dcds_BDDb_im", -100., 100.);
        AddSU3BreakingParameter("delta_G3t_cds_BDDb_re", "G3t_css_BDDb_re");
        AddSU3BreakingParameter("delta_G3t_cds_BDDb_im", "G3t_css_BDDb_im");
    }
    else if (channel == "Bddpsdms")
    {
        // b → c(c̄d), spectator d
        std::vector<std::string> params = {
            "delta_A2t_cscd_BDbD_re", "delta_A2t_cscd_BDbD_im",
            "delta_G3t_csd_BDDb_re", "delta_G3t_csd_BDDb_im"};
        channelParameters[channel] = params;

        AddSU3BreakingParameter("delta_A2t_cscd_BDbD_re", "A2t_cscs_BDbD_re");
        AddSU3BreakingParameter("delta_A2t_cscd_BDbD_im", "A2t_cscs_BDbD_im");
        AddSU3BreakingParameter("delta_G3t_csd_BDDb_re", "G3t_css_BDDb_re");
        AddSU3BreakingParameter("delta_G3t_csd_BDDb_im", "G3t_css_BDDb_im");
    }
    else if (channel == "Bddpsdm")
    {
        // b → c(c̄s), spectator d
        std::vector<std::string> params = {
            "delta_E1t_sccd_BDDb_re", "delta_E1t_sccd_BDDb_im",
            "delta_G1t_scd_BDDb_re", "delta_G1t_scd_BDDb_im"};
        channelParameters[channel] = params;

        AddSU3BreakingParameter("delta_E1t_sccd_BDDb_re", "E1t_sccs_BDDb_re");
        AddSU3BreakingParameter("delta_E1t_sccd_BDDb_im", "E1t_sccs_BDDb_im");
        AddSU3BreakingParameter("delta_G1t_scd_BDDb_re", "G1t_scs_BDDb_re");
        AddSU3BreakingParameter("delta_G1t_scd_BDDb_im", "G1t_scs_BDDb_im");
    }
    else if (channel == "Bddpdm")
    {
        // b → c(c̄d), spectator d
        std::vector<std::string> params = {
            "delta_E1t_dccs_BDDb_re", "delta_E1t_dccs_BDDb_im",
            "delta_E1t_dccd_BDDb_re", "delta_E1t_dccd_BDDb_im",
            "delta_A2t_cdcs_BDbD_re", "delta_A2t_cdcs_BDbD_im",
            "delta_A2t_cdcd_BDbD_re", "delta_A2t_cdcd_BDbD_im",
            "delta_G1t_dcs_BDDb_re", "delta_G1t_dcs_BDDb_im",
            "delta_G3t_cds_BDDb_re", "delta_G3t_cds_BDDb_im",
            "delta_G1t_dcd_BDDb_re", "delta_G1t_dcd_BDDb_im",
            "delta_G3t_cdd_BDDb_re", "delta_G3t_cdd_BDDb_im"};
        channelParameters[channel] = params;

        addSU3BreakingParameter("delta_E1t_dccs_BDDb_re", "E1t_sccs_BDDb_re");
        addSU3BreakingParameter("delta_E1t_dccs_BDDb_im", "E1t_sccs_BDDb_im");
        addSU3BreakingParameter("delta_E1t_dccd_BDDb_re", "E1t_dccs_BDDb_re");
        addSU3BreakingParameter("delta_E1t_dccd_BDDb_im", "E1t_dccs_BDDb_im");
        addSU3BreakingParameter("delta_A2t_cdcs_BDbD_re", "A2t_cscs_BDbD_re");
        addSU3BreakingParameter("delta_A2t_cdcs_BDbD_im", "A2t_cscs_BDbD_im");
        addSU3BreakingParameter("delta_A2t_cdcd_BDbD_re", "A2t_cdcs_BDbD_re");
        addSU3BreakingParameter("delta_A2t_cdcd_BDbD_im", "A2t_cdcs_BDbD_im");
        addSU3BreakingParameter("delta_G1t_dcs_BDDb_re", "G1t_scs_BDDb_re");
        addSU3BreakingParameter("delta_G1t_dcs_BDDb_im", "G1t_scs_BDDb_im");
        addSU3BreakingParameter("delta_G3t_cds_BDDb_re", "G3t_css_BDDb_re");
        addSU3BreakingParameter("delta_G3t_cds_BDDb_im", "G3t_css_BDDb_im");
        addSU3BreakingParameter("delta_G1t_dcd_BDDb_re", "G1t_dcs_BDDb_re");
        addSU3BreakingParameter("delta_G1t_dcd_BDDb_im", "G1t_dcs_BDDb_im");
        addSU3BreakingParameter("delta_G3t_cdd_BDDb_re", "G3t_cds_BDDb_re");
        addSU3BreakingParameter("delta_G3t_cdd_BDDb_im", "G3t_cds_BDDb_im");
    }
    else if (channel == "Bdd0d0b")
    {
        // b → c(c̄d), spectator d
        std::vector<std::string> params = {
            "delta_A2t_cdcs_BDbD_re", "delta_A2t_cdcs_BDbD_im",
            "delta_A2t_cdcd_BDbD_re", "delta_A2t_cdcd_BDbD_im",
            "delta_dP3EW_ucd_BDbD_re", "delta_dP3EW_ucd_BDbD_im",
            "delta_A2_dcdd_BDDb_re", "delta_A2_dcdd_BDDb_im",
            "delta_G3t_cds_BDDb_re", "delta_G3t_cds_BDDb_im",
            "delta_G3t_cdd_BDDb_re", "delta_G3t_cdd_BDDb_im"};
        channelParameters[channel] = params;

        AddSU3BreakingParameter("delta_A2t_cdcs_BDbD_re", "A2t_cscs_BDbD_re");
        AddSU3BreakingParameter("delta_A2t_cdcs_BDbD_im", "A2t_cscs_BDbD_im");
        AddSU3BreakingParameter("delta_A2t_cdcd_BDbD_re", "A2t_cdcs_BDbD_re");
        AddSU3BreakingParameter("delta_A2t_cdcd_BDbD_im", "A2t_cdcs_BDbD_im");
        addSU3BreakingParameter("delta_dP3EW_ucd_BDbD_re", "dP3EW_ucs_BDbD_re");
        addSU3BreakingParameter("delta_dP3EW_ucd_BDbD_im", "dP3EW_ucs_BDbD_im");
        addSU3BreakingParameter("delta_A2_dcdd_BDDb_re", "A2_dcds_BDDb_re");
        addSU3BreakingParameter("delta_A2_dcdd_BDDb_im", "A2_dcds_BDDb_im");
        AddSU3BreakingParameter("delta_G3t_cds_BDDb_re", "G3t_css_BDDb_re");
        AddSU3BreakingParameter("delta_G3t_cds_BDDb_im", "G3t_css_BDDb_im");
        AddSU3BreakingParameter("delta_G3t_cdd_BDDb_re", "G3t_cds_BDDb_re");
        AddSU3BreakingParameter("delta_G3t_cdd_BDDb_im", "G3t_cds_BDDb_im");
    }
    else if (channel == "Bpdpd0b")
    {
        // b → c(c̄d), spectator u
        std::vector<std::string> params = {
            "delta_E1t_dccs_BDDb_re", "delta_E1tt_dccs_BDDb_im",
            "delta_E1t_dccd_BDDb_re", "delta_E1t_dccd_BDDb_im",
            "dP1EW_dcu_BDDb_re", "dP1EW_dcu_BDDb_im",
            "A1_dcdd_BDDb_re", "A1_dcdd_BDDb_im",
            "delta_G1t_dcd_BDDb_re", "delta_G1t_dcd_BDDb_im"};
        channelParameters[channel] = params;

        addSU3BreakingParameter("delta_E1t_dccs_BDDb_re", "E1t_sccs_BDDb_re");
        addSU3BreakingParameter("delta_E1t_dccs_BDDb_im", "E1t_sccs_BDDb_im");
        addSU3BreakingParameter("delta_E1t_dccd_BDDb_re", "E1t_dccs_BDDb_re");
        addSU3BreakingParameter("delta_E1t_dccd_BDDb_im", "E1t_dccs_BDDb_im");
        addAmplitudeParameter("dP1EW_dcu_BDDb_re", -ewp_limit, ewp_limit);
        addAmplitudeParameter("dP1EW_dcu_BDDb_im", -ewp_limit, ewp_limit);
        addAmplitudeParameter("A1_dcdd_BDDb_re", -100., 100.);
        addAmplitudeParameter("A1_dcdd_BDDb_im", -100., 100.);
        addSU3BreakingParameter("delta_G1t_dcd_BDDb_re", "G1t_dcs_BDDb_re");
        addSU3BreakingParameter("delta_G1t_dcd_BDDb_im", "G1t_dcs_BDDb_im");
    }
    else if (channel == "Bpdspd0b")
    {
        // b → c(c̄s), spectator u
        std::vector<std::string> params = {
            "delta_E1t_sccd_BDDb_re", "delta_E1t_sccd_BDDb_im",
            "delta_dP1EW_scu_BDDb_re", "delta_dP1EW_scu_BDDb_im",
            "delta_A1_scdd_BDDb_re", "delta_A1_scdd_BDDb_im",
            "delta_G1t_scd_BDDb_re", "delta_G1t_scd_BDDb_im"};
        channelParameters[channel] = params;

        addSU3BreakingParameter("delta_E1t_sccd_BDDb_re", "E1t_sccs_BDDb_re");
        addSU3BreakingParameter("delta_E1t_sccd_BDDb_im", "E1t_sccs_BDDb_im");
        addSU3BreakingParameter("delta_dP1EW_scu_BDDb_re", "dP1EW_dcu_BDDb_re");
        addSU3BreakingParameter("delta_dP1EW_scu_BDDb_im", "dP1EW_dcu_BDDb_im");
        addSU3BreakingParameter("delta_A1_scdd_BDDb_re", "A1_dcdd_BDDb_re");
        addSU3BreakingParameter("delta_A1_scdd_BDDb_im", "A1_dcdd_BDDb_im");
        addSU3BreakingParameter("delta_G1t_scd_BDDb_re", "G1t_scs_BDDb_re");
        addSU3BreakingParameter("delta_G1t_scd_BDDb_im", "G1t_scs_BDDb_im");
    }
    else
    {
        // Handle the case where the channel is not recognized
        std::cerr << "Error: Unrecognized channel \"" << channel << "\" in DefineParameters." << std::endl;
        throw std::runtime_error("Unrecognized channel: " + channel);
    }

    std::cout << "Number of parameters: " << GetNParameters() << std::endl;
    SetPriorConstantAll();
}

std::map<std::string, double> parameterValues;
// function that given the channels and the string inside the map channelParameters adds every parameter name to a map <string, double> and returns said map.
std::map<std::string, double> goldenmodesB::DeclareParameters()
{

    // Ensure channelParameters is populated
    for (const auto &channel : channelNamesSU3)
    {
        // Call DefineParameters once per channel
        if (channelParameters.find(channel) == channelParameters.end())
        {
            goldenmodesB::DefineParameters(channel);
        }

        // Add parameters to parameterValues
        auto it = channelParameters.find(channel);
        if (it != channelParameters.end())
        {
            for (const auto &param : it->second)
            {
                parameterValues[param] = 1.; // initialize
                cout << "parameter " << param << " added with value: " << parameterValues[param] << endl;
            }
        }
        else
        {
            std::cerr << "Channel " << channel << " not found in channelParameters." << std::endl;
        }
    }

    cout << "Declare Parameters called correctly" << endl;
    return parameterValues;
}

TComplex goldenmodesB::getPar(const std::string &baseName) const
{
    // Look for the real and imaginary parts in the parameterValues map
    auto it_real = parameterValues.find(baseName + "_re");
    auto it_imag = parameterValues.find(baseName + "_im");

    if (it_real != parameterValues.end() && it_imag != parameterValues.end())
    {
        // Construct a complex number from the real and imaginary parts
        return TComplex(it_real->second, it_imag->second);
    }
    else
    {
        throw std::runtime_error("Error: Real or imaginary part for parameter " + baseName + " not found.");
    }
}

// Setter function: sets the value for a given parameter in the map
void goldenmodesB::SetParameterValue(const std::string &paramName, double value)
{
    parameterValues[paramName] = value; // Insert or update the value for the given parameter name
}

//---------------------------------------------------------------------

double goldenmodesB::getParameterValue(const std::string &paramName) const
{
    // Find the parameter in the parameterValues map
    auto it = parameterValues.find(paramName);

    // If the parameter is found, return its value
    if (it != parameterValues.end())
    {
        return it->second;
    }
    else
    {
        // If the parameter is not found, throw an error or return a default value
        std::cerr << "Error: Parameter " << paramName << " not found in parameterValues map." << std::endl;
        throw std::runtime_error("Parameter not found: " + paramName);
    }
}

//----------------------------------------------------------

// compute decay amplitude for each channel
void goldenmodesB::compute_decay_amplitudes(const std::string &channel, bool conjugate)
{
    // cout << "computing decay amplitude for channel " << channel << endl;
    amplitudes[channel] = TComplex(0.0, 0.0);

    Parameter amp;
    Parameter amp_0;
    Parameter amp_paral;
    Parameter amp_perp;
    // Get the CKM elements for the current channel, apply conjugation based on the 'conjugate' flag
    TComplex lam_bs_c = conjugate ? ckm.getVcb() * TComplex::Conjugate(ckm.getVcs()) : ckm.getVcs() * TComplex::Conjugate(ckm.getVcb());
    TComplex lam_bs_u = conjugate ? ckm.getVub() * TComplex::Conjugate(ckm.getVus()) : ckm.getVus() * TComplex::Conjugate(ckm.getVub());
    TComplex lam_bd_c = conjugate ? ckm.getVcb() * TComplex::Conjugate(ckm.getVcd()) : ckm.getVcd() * TComplex::Conjugate(ckm.getVcb());
    TComplex lam_bd_u = conjugate ? ckm.getVub() * TComplex::Conjugate(ckm.getVud()) : ckm.getVud() * TComplex::Conjugate(ckm.getVub());

    if (channel == "Bdjpsik0")
    {
        // Bd→J/ψ K⁰: b→c(c̄s), spectator d
        amp = lam_bs_c * getPar("E2t_ccsd_BJPSIP") - lam_bs_u * getPar("G2t_scd_BJPSIP");
        amplitudes[channel] = amp;
    }
    else if (channel == "Bdjpsip0")
    {
        // Bd→J/ψ π⁰: b→c(c̄d), spectator d
        amp = (lam_bd_c * (getPar("E2t_ccdd_BJPSIP") + getPar("dP4EW_ucd_BdPJPSI")) -
               lam_bd_u * (getPar("EA2_ddcd_BPJPSI") + getPar("dP4EW_ucd_BdPJPSI") + getPar("G2t_dcd_BJPSIP"))) /
              sqrt(2.);
        amplitudes[channel] = amp;
    }
    else if (channel == "Bdjpsiet8")
    {
        // Bd→J/ψ eta_8: b→c(c̄d), spectator d
        amp = ((lam_bd_c * (getPar("E2t_ccdd_BJPSIP") + getPar("dP4EW_ucd_BPJPSI") + 2. * getPar("EA2t_ccdd_BJPSIP") - 2. * getPar("EA2t_ccsd_BJPSIP"))) +
               lam_bd_u * (getPar("EA2_ddcd_BPJPSI") + getPar("dP4EW_ucd_BPJPSI") - getPar("G2t_dcd_BJPSIP") - 2. * getPar("G4t_cdd_BJPSIP") + 2. * getPar("G4t_csd_BJPSIP"))) /
              sqrt(6.);
        amplitudes[channel] = amp;
    }
    else if (channel == "Bdjpsiet1")
    {
        // Bd→J/ψ eta_1: b→c(c̄d), spectator d
        amp = ((lam_bd_c * (getPar("E2t_ccdd_BJPSIP") + getPar("dP4EW_ucd_BPJPSI") + 2. * getPar("EA2t_ccdd_BJPSIP") + getPar("EA2t_ccsd_BJPSIP"))) +
               lam_bd_u * (getPar("EA2_ddcd_BPJPSI") + getPar("dP4EW_ucd_BPJPSI") - getPar("G2t_dcd_BJPSIP") - 2. * getPar("G4t_cdd_BJPSIP") - getPar("G4t_csd_BJPSIP"))) /
              sqrt(3.);
        amplitudes[channel] = amp;
    }
    else if (channel == "Bpjpsikp")
    {
        // B⁺→J/ψ K⁺: b→c(c̄s), spectator u
        amp = lam_bs_c * (getPar("E2t_ccsd_BJPSIP") + getPar("dP2EW_scu_Bpjpsikp")) +
              lam_bs_u * (getPar("EA1_sdcd_BPJPSI") + getPar("dP2EW_scu_Bpjpsikp") - getPar("G2t_scd_BJPSIP"));
        amplitudes[channel] = amp;
    }
    else if (channel == "Bpjpsipp")
    {
        // B⁺→J/ψ π⁺: b→c(c̄d), spectator u
        amp = lam_bd_c * (getPar("E2t_ccdd_BJPSIP") + getPar("dP2EW_dcu_Bpjpsipp")) +
              lam_bd_u * (getPar("EA1_ddcd_BPJPSI") + getPar("dP2EW_dcu_Bpjpsipp") - getPar("G2t_dcd_BJPSIP"));
        amplitudes[channel] = amp;
    }
    else if (channel == "Bsjpsip0")
    {
        // Bs→J/ψ π⁰: b→c(c̄s), spectator s
        amp = -(lam_bs_c * getPar("dP4EW_ucs_BPJPSI") + lam_bs_u * (getPar("EA2_ddcs_BPJPSI") + getPar("dP4EW_ucs_BPJPSI"))) / sqrt(2.);
        amplitudes[channel] = amp;
    }
    else if (channel == "Bsjpsik0b")
    {
        // Bs→J/ψ \bar{K}⁰: b→c(c̄d), spectator s
        amp = lam_bd_c * getPar("E2t_ccds_BJPSIP") - lam_bd_u * getPar("G2t_dcs_BJPSIP");
        amplitudes[channel] = amp;
    }
    else if (channel == "Bsjpsiet8")
    {
        // Bs→J/ψ eta_8: b→c(c̄s), spectator s
        amp = ((lam_bs_c * (-2. * getPar("E2t_ccss_BJPSIP") + getPar("dP4EW_ucs_BPJPSI") + 2. * getPar("EA2t_ccds_BJPSIP") - 2. * getPar("EA2t_ccss_BJPSIP"))) +
               lam_bs_u * (getPar("EA2_ddcs_BPJPSI") + getPar("dP4EW_ucs_BPJPSI") + 2. * getPar("G2t_scs_BJPSIP") - 2. * getPar("G4t_cds_BJPSIP") + 2. * getPar("G4t_css_BJPSIP"))) /
              sqrt(6.);
        amplitudes[channel] = amp;
    }
    else if (channel == "Bsjpsiet1")
    {
        // Bs→J/ψ eta_1: b→c(c̄s), spectator s
        amp = ((lam_bs_c * (getPar("E2t_ccss_BJPSIP") + getPar("dP4EW_ucs_BPJPSI") + 2. * getPar("EA2t_ccds_BJPSIP") + getPar("EA2t_ccss_BJPSIP"))) +
               lam_bs_u * (getPar("EA2_ddcs_BPJPSI") + getPar("dP4EW_ucs_BPJPSI") - getPar("G2t_scs_BJPSIP") - 2. * getPar("G4t_cds_BJPSIP") - getPar("G4t_css_BJPSIP"))) /
              sqrt(3.);
        amplitudes[channel] = amp;
    }
    else if (channel == "Bsjpsiph")
    {
        // Bs→J/ψ φ: b→c(c̄s), spectator s
        amp_0 = -lam_bs_c * (getPar("E2t_ccss_BJPSIV_0") + getPar("EA2t_ccss_BJPSIV_0")) +
                lam_bs_u * (getPar("G2t_scs_BJPSIV_0") + getPar("G4t_css_BJPSIV_0"));
        amp_paral = -lam_bs_c * (getPar("E2t_ccss_BJPSIV_paral") + getPar("EA2t_ccss_BJPSIV_paral")) +
                    lam_bs_u * (getPar("G2t_scs_BJPSIV_paral") + getPar("G4t_css_BJPSIV_paral"));
        amp_perp = -lam_bs_c * (getPar("E2t_ccss_BJPSIV_perp") + getPar("EA2t_ccss_BJPSIV_perp")) +
                   lam_bs_u * (getPar("G2t_scs_BJPSIV_perp") + getPar("G4t_css_BJPSIV_perp"));
        amplitudes[channel + "_0"] = amp_0;
        amplitudes[channel + "_paral"] = amp_paral;
        amplitudes[channel + "_perp"] = amp_perp;
    }
    else if (channel == "Bsjpsiom")
    {
        // Bs→J/ψ ω: b→c(c̄s), spectator s
        amp_0 = (lam_bs_c * (2. * getPar("EA2t_ccds_BJPSIV_0") + getPar("dP4EW_ucs_BVJPSI_0")) +
                 lam_bs_u * (getPar("EA2_ddcs_BVJPSI_0") + getPar("dP4EW_ucs_BVJPSI_0") - 2. * getPar("G4t_cds_BJPSIV_0"))) /
                sqrt(2.);
        amp_paral = (lam_bs_c * (2. * getPar("EA2t_ccds_BJPSIV_paral") + getPar("dP4EW_ucs_BVJPSI_paral")) +
                     lam_bs_u * (getPar("EA2_ddcs_BVJPSI_paral") + getPar("dP4EW_ucs_BVJPSI_paral") - 2. * getPar("G4t_cds_BJPSIV_paral"))) /
                    sqrt(2.);
        amp_perp = (lam_bs_c * (2. * getPar("EA2t_ccds_BJPSIV_perp") + getPar("dP4EW_ucs_BVJPSI_perp")) +
                    lam_bs_u * (getPar("EA2_ddcs_BVJPSI_perp") + getPar("dP4EW_ucs_BVJPSI_perp") - 2. * getPar("G4t_cds_BJPSIV_perp"))) /
                   sqrt(2.);
        amplitudes[channel + "_0"] = amp_0;
        amplitudes[channel + "_paral"] = amp_paral;
        amplitudes[channel + "_perp"] = amp_perp;
    }
    else if (channel == "Bsjpsikbst")
    {
        // Bs→J/ψ \bar{K}*: b→c(c̄d), spectator s
        amp_0 = -lam_bd_c * getPar("E2t_ccds_BJPSIV_0") - lam_bd_u * getPar("G2t_dcs_BJPSIV_0");
        amp_paral = -lam_bd_c * getPar("E2t_ccds_BJPSIV_paral") - lam_bd_u * getPar("G2t_dcs_BJPSIV_paral");
        amp_perp = -lam_bd_c * getPar("E2t_ccds_BJPSIV_perp") - lam_bd_u * getPar("G2t_dcs_BJPSIV_perp");
        amplitudes[channel + "_0"] = amp_0;
        amplitudes[channel + "_paral"] = amp_paral;
        amplitudes[channel + "_perp"] = amp_perp;
    }
    else if (channel == "Bsjpsirho0")
    {
        // Bs→J/ψ \rho⁰: b→c(c̄s), spectator s
        amp_0 = -(lam_bs_c * getPar("dP4EW_ucs_BVJPSI_0") + lam_bs_u * (getPar("EA2_ddcs_BVJPSI_0") + getPar("dP4EW_ucs_BVJPSI_0"))) / sqrt(2.);
        amp_paral = -(lam_bs_c * getPar("dP4EW_ucs_BVJPSI_paral") + lam_bs_u * (getPar("EA2_ddcs_BVJPSI_paral") + getPar("dP4EW_ucs_BVJPSI_paral"))) / sqrt(2.);
        amp_perp = -(lam_bs_c * getPar("dP4EW_ucs_BVJPSI_perp") + lam_bs_u * (getPar("EA2_ddcs_BVJPSI_perp") + getPar("dP4EW_ucs_BVJPSI_perp"))) / sqrt(2.);
        amplitudes[channel + "_0"] = amp_0;
        amplitudes[channel + "_paral"] = amp_paral;
        amplitudes[channel + "_perp"] = amp_perp;
    }
    else if (channel == "Bdjpsiom")
    {
        // Bd→J/ψ ω: b→c(c̄d), spectator d
        amp_0 = (lam_bd_c * (getPar("E2t_ccdd_BJPSIV_0") + 2. * getPar("EA2t_ccdd_BJPSIV_0") + getPar("dP4EW_ucd_BVJPSI_0")) +
                 lam_bd_u * (getPar("EA2_ddcd_BVJPSI_0") + getPar("dP4EW_ucd_BVJPSI_0") - getPar("G2t_dcd_BJPSIV_0") - 2. * getPar("G4t_cdd_BJPSIV_0"))) /
                sqrt(2.);
        amp_paral = (lam_bd_c * (getPar("E2t_ccdd_BJPSIV_paral") + 2. * getPar("EA2t_ccdd_BJPSIV_paral") + getPar("dP4EW_ucd_BVJPSI_paral")) +
                     lam_bd_u * (getPar("EA2_ddcd_BVJPSI_paral") + getPar("dP4EW_ucd_BVJPSI_paral") - getPar("G2t_dcd_BJPSIV_paral") - 2. * getPar("G4t_cdd_BJPSIV_paral"))) /
                    sqrt(2.);
        amp_perp = (lam_bd_c * (getPar("E2t_ccdd_BJPSIV_perp") + 2. * getPar("EA2t_ccdd_BJPSIV_perp") + getPar("dP4EW_ucd_BVJPSI_perp")) +
                    lam_bd_u * (getPar("EA2_ddcd_BVJPSI_perp") + getPar("dP4EW_ucd_BVJPSI_perp") - getPar("G2t_dcd_BJPSIV_perp") - 2. * getPar("G4t_cdd_BJPSIV_perp"))) /
                   sqrt(2.);
        amplitudes[channel + "_0"] = amp_0;
        amplitudes[channel + "_paral"] = amp_paral;
        amplitudes[channel + "_perp"] = amp_perp;
    }
    else if (channel == "Bdjpsikst")
    {
        // Bd→J/ψ K*: b→c(c̄s), spectator d
        amp_0 = lam_bs_c * getPar("E2t_ccsd_BJPSIV_0") - lam_bs_u * getPar("G2t_scd_BJPSIV_0");
        amp_paral = lam_bs_c * getPar("E2t_ccsd_BJPSIV_paral") - lam_bs_u * getPar("G2t_scd_BJPSIV_paral");
        amp_perp = lam_bs_c * getPar("E2t_ccsd_BJPSIV_perp") - lam_bs_u * getPar("G2t_scd_BJPSIV_perp");
        amplitudes[channel + "_0"] = amp_0;
        amplitudes[channel + "_paral"] = amp_paral;
        amplitudes[channel + "_perp"] = amp_perp;
    }
    else if (channel == "Bdjpsirho0")
    {
        // Bd→J/ψ ρ: b→c(c̄d), spectator d
        amp_0 = (lam_bd_c * (getPar("E2t_ccdd_BJPSIV_0") + getPar("dP4EW_ucd_BVJPSI_0")) -
                 lam_bd_u * (getPar("EA2_ddcd_BVJPSI_0") + getPar("dP4EW_ucd_BVJPSI_0") + getPar("G2t_dcd_BJPSIV_0"))) /
                sqrt(2.);
        amp_paral = (lam_bd_c * (getPar("E2t_ccdd_BJPSIV_paral") + getPar("dP4EW_ucd_BVJPSI_paral")) -
                     lam_bd_u * (getPar("EA2_ddcd_BVJPSI_paral") + getPar("dP4EW_ucd_BVJPSI_paral") + getPar("G2t_dcd_BJPSIV_paral"))) /
                    sqrt(2.);
        amp_perp = (lam_bd_c * (getPar("E2t_ccdd_BJPSIV_perp") + getPar("dP4EW_ucd_BVJPSI_perp")) -
                    lam_bd_u * (getPar("EA2_ddcd_BVJPSI_perp") + getPar("dP4EW_ucd_BVJPSI_perp") + getPar("G2t_dcd_BJPSIV_perp"))) /
                   sqrt(2.);
        amplitudes[channel + "_0"] = amp_0;
        amplitudes[channel + "_paral"] = amp_paral;
        amplitudes[channel + "_perp"] = amp_perp;
    }
    else if (channel == "Bdjpsiph")
    {
        // Bd→J/ψ \phi: b→c(c̄d), spectator d
        amp_0 = -lam_bd_c * getPar("EA2t_ccsd_BJPSIV_0") - lam_bd_u * getPar("G4t_csd_BJPSIV_0");
        amp_paral = -lam_bd_c * getPar("EA2t_ccsd_BJPSIV_paral") - lam_bd_u * getPar("G4t_csd_BJPSIV_paral");
        amp_perp = -lam_bd_c * getPar("EA2t_ccsd_BJPSIV_perp") - lam_bd_u * getPar("G4t_csd_BJPSIV_perp");
        amplitudes[channel + "_0"] = amp_0;
        amplitudes[channel + "_paral"] = amp_paral;
        amplitudes[channel + "_perp"] = amp_perp;
    }
    else if (channel == "Bpjpsirhop")
    {
        // B⁺→J/ψ ρ⁺: b→c(c̄d), spectator u
        amp_0 = lam_bd_c * (getPar("E2t_ccdd_BJPSIV_0") + getPar("dP2EW_dcu_BJPSIV_0")) +
                lam_bd_u * (getPar("EA1_ddcd_BVJPSI_0") + getPar("dP2EW_dcu_BJPSIV_0") - getPar("G2t_dcd_BJPSIV_0"));
        amp_paral = lam_bd_c * (getPar("E2t_ccdd_BJPSIV_paral") + getPar("dP2EW_dcu_BJPSIV_paral")) +
                    lam_bd_u * (getPar("EA1_ddcd_BVJPSI_paral") + getPar("dP2EW_dcu_BJPSIV_paral") - getPar("G2t_dcd_BJPSIV_paral"));
        amp_perp = lam_bd_c * (getPar("E2t_ccdd_BJPSIV_perp") + getPar("dP2EW_dcu_BJPSIV_perp")) +
                   lam_bd_u * (getPar("EA1_ddcd_BVJPSI_perp") + getPar("dP2EW_dcu_BJPSIV_perp") - getPar("G2t_dcd_BJPSIV_perp"));
        amplitudes[channel + "_0"] = amp_0;
        amplitudes[channel + "_paral"] = amp_paral;
        amplitudes[channel + "_perp"] = amp_perp;
    }
    else if (channel == "Bpjpsikstp")
    {
        // B⁺→J/ψ K*⁺: b→c(c̄s), spectator u
        amp_0 = lam_bs_c * (getPar("E2t_ccsd_BJPSIV_0") + getPar("dP2EW_scu_BJPSIV_0")) +
                lam_bs_u * (getPar("EA1_sdcd_BVJPSI_0") + getPar("dP2EW_scu_BJPSIV_0") - getPar("G2t_scd_BJPSIV_0"));
        amp_paral = lam_bs_c * (getPar("E2t_ccsd_BJPSIV_paral") + getPar("dP2EW_scu_BJPSIV_paral")) +
                    lam_bs_u * (getPar("EA1_sdcd_BVJPSI_paral") + getPar("dP2EW_scu_BJPSIV_paral") - getPar("G2t_scd_BJPSIV_paral"));
        amp_perp = lam_bs_c * (getPar("E2t_ccsd_BJPSIV_perp") + getPar("dP2EW_scu_BJPSIV_perp")) +
                   lam_bs_u * (getPar("EA1_sdcd_BVJPSI_perp") + getPar("dP2EW_scu_BJPSIV_perp") - getPar("G2t_scd_BJPSIV_perp"));
        amplitudes[channel + "_0"] = amp_0;
        amplitudes[channel + "_paral"] = amp_paral;
        amplitudes[channel + "_perp"] = amp_perp;
    }
    else if (channel == "Bsdpsdms")
    {
        // b → c(c̄s), spectator s
        amp = lam_bs_c * (getPar("E1t_sccs_BDDb") + getPar("A2t_cscs_BDbD")) - lam_bs_u * (getPar("G1t_scs_BDDb") + getPar("G3t_css_BDDb"));
        amplitudes[channel] = amp;
    }
    else if (channel == "Bsdpdms")
    {
        // b → c(c̄d), spectator s
        amp = lam_bd_c * (getPar("E1t_dccs_BDDb")) - lam_bs_u * (getPar("G1t_dcs_BDDb"));
        amplitudes[channel] = amp;
    }
    else if (channel == "Bsdpdm")
    {
        // b → c(c̄s), spectator s
        amp = lam_bs_c * (getPar("A2t_cdcs_BDbD")) - lam_bs_u * (getPar("G3t_cds_BDDb"));
        amplitudes[channel] = amp;
    }
    else if (channel == "Bsd0d0b")
    {
        // b → c(c̄s), spectator s
        amp = -lam_bs_c * (getPar("A2t_cdcs_BDbD") + getPar("dP3EW_ucs_BDbD")) -
              lam_bs_u * (getPar("A2_dcds_BDDb") + getPar("dP3EW_ucs_BDbD") - getPar("G3t_cds_BDDb"));
        amplitudes[channel] = amp;
    }
    else if (channel == "Bddpsdms")
    {
        // b → c(c̄d), spectator d
        amp = lam_bd_c * (getPar("A2t_cscd_BDbD")) - lam_bd_u * (getPar("G3t_csd_BDDb"));
        amplitudes[channel] = amp;
    }
    else if (channel == "Bddpsdm")
    {
        // b → c(c̄s), spectator d
        amp = lam_bs_c * (getPar("E1t_sccd_BDDb")) - lam_bs_u * (getPar("G1t_scd_BDDb"));
        amplitudes[channel] = amp;
    }
    else if (channel == "Bddpdm")
    {
        // b → c(c̄d), spectator d
        amp = lam_bd_c * (getPar("E1t_dccd_BDDb") + getPar("A2t_cdcd_BDbD")) -
              lam_bd_u * (getPar("G1t_dcd_BDDb") + getPar("G3t_cdd_BDDb"));
        amplitudes[channel] = amp;
    }
    else if (channel == "Bdd0d0b")
    {
        // b → c(c̄d), spectator d
        amp = -lam_bd_c * (getPar("A2t_cdcd_BDbD") + getPar("dP3EW_ucd_BDbD")) -
              lam_bd_u * (getPar("A2_dcdd_BDDb") + getPar("dP3EW_ucd_BDbD") - getPar("G3t_cdd_BDDb"));
        amplitudes[channel] = amp;
    }
    else if (channel == "Bpdpd0b")
    {
        // b → c(c̄d), spectator u
        amp = +lam_bd_c * (getPar("E1t_dccd_BDDb") + getPar("dP1EW_dcu_BDDb")) +
              lam_bd_u * (getPar("A1_dcdd_BDDb") + getPar("dP1EW_dcu_BDDb") - getPar("G1t_dcd_BDDb"));
        amplitudes[channel] = amp;
    }
    else if (channel == "Bpdspd0b")
    {
        // b → c(c̄s), spectator u
        amp = +lam_bs_c * (getPar("E1t_sccd_BDDb") + getPar("dP1EW_scu_BDDb")) +
              lam_bs_u * (getPar("A1_scdd_BDDb") + getPar("dP1EW_scu_BDDb") - getPar("G1t_scd_BDDb"));
        amplitudes[channel] = amp;
    }
    else
    {
        cout << "WARNING: amplitude for channel " << channel << " not found" << endl;
    }
}

// Getter for amplitudes
Parameter goldenmodesB::get_amplitude(const std::string &channel)
{
    // First, check if the amplitude is already computed
    auto it = amplitudes.find(channel);
    if (it != amplitudes.end())
    {
        return it->second;
    }

    // Extract the base channel name (if it's a polarized component)
    std::string base_channel = channel;
    if (channel.find("_0") != std::string::npos)
    {
        base_channel = channel.substr(0, channel.find("_0"));
    }
    else if (channel.find("_paral") != std::string::npos)
    {
        base_channel = channel.substr(0, channel.find("_paral"));
    }
    else if (channel.find("_perp") != std::string::npos)
    {
        base_channel = channel.substr(0, channel.find("_perp"));
    }

    // Compute the decay amplitude for the base channel
    compute_decay_amplitudes(base_channel, false);

    // Check if the specific requested amplitude exists
    it = amplitudes.find(channel);
    if (it == amplitudes.end())
    {
        std::cerr << "ERROR: Amplitude not found for channel: " << channel << std::endl;
        std::cerr << "Available amplitudes: ";
        for (const auto &pair : amplitudes)
        {
            std::cerr << pair.first << " ";
        }
        std::cerr << std::endl;
        throw std::runtime_error("Amplitude not found for channel: " + channel);
    }

    return it->second;
}

// Getter for conjugated amplitudes
Parameter goldenmodesB::get_conjugate_amplitude(const std::string &channel)
{
    // First, check if the conjugate amplitude already exists
    auto it = amplitudes.find(channel);
    if (it != amplitudes.end())
    {
        return it->second;
    }

    // Extract the base channel name (if it's a polarized component)
    std::string base_channel = channel;
    if (channel.find("_0") != std::string::npos)
    {
        base_channel = channel.substr(0, channel.find("_0"));
    }
    else if (channel.find("_paral") != std::string::npos)
    {
        base_channel = channel.substr(0, channel.find("_paral"));
    }
    else if (channel.find("_perp") != std::string::npos)
    {
        base_channel = channel.substr(0, channel.find("_perp"));
    }

    // Compute the conjugate decay amplitude for the base channel
    compute_decay_amplitudes(base_channel, true);
    // Check if the specific conjugate amplitude exists
    it = amplitudes.find(channel);
    if (it == amplitudes.end())
    {
        std::cerr << "ERROR: Conjugate amplitude not found for channel: " << channel << std::endl;
        std::cerr << "Available amplitudes: ";
        for (const auto &pair : amplitudes)
        {
            std::cerr << pair.first << " ";
        }
        std::cerr << std::endl;
        throw std::runtime_error("Conjugate amplitude not found for channel: " + channel);
    }

    return it->second;
}

// Destructor
goldenmodesB::~goldenmodesB()
{
    // Ensure any dynamically allocated resources are properly cleaned
    // delete histos;  // Uncomment if histos is dynamically allocated
}

//-----------------------------------------------------------
// Helper function to parse a decay channel name
std::pair<std::string, std::pair<std::string, std::string>>
goldenmodesB::parseChannel(const std::string &channel) const
{
    // Identify the B meson part
    std::string bMeson;
    if (channel.rfind("Bp", 0) == 0)
    {
        bMeson = "Bp";
    }
    else if (channel.rfind("Bd", 0) == 0)
    {
        bMeson = "Bd";
    }
    else if (channel.rfind("Bs", 0) == 0)
    {
        bMeson = "Bs";
    }
    else
    {
        throw std::runtime_error("Error in parseChannel: Unknown B meson in channel: " + channel);
    }

    std::string remaining = channel.substr(bMeson.length());

    // Check if it's a J/psi channel or a DD channel
    const std::string jpsiIdentifier = "jpsi";
    if (remaining.rfind(jpsiIdentifier, 0) == 0)
    {
        // J/psi channel: B→J/ψ X
        // Extract the second final-state meson
        remaining = remaining.substr(jpsiIdentifier.length());
        std::string meson2;

        // Look for a valid meson name in mesonMasses
        for (const auto &meson : mesonMasses)
        {
            if (remaining.rfind(meson.first, 0) == 0)
            {
                meson2 = meson.first;
                break;
            }
        }

        if (meson2.empty())
        {
            throw std::runtime_error("Error in parseChannel: Unable to parse second final-state meson in channel: " + channel);
        }

        return {bMeson, {jpsiIdentifier, meson2}};
    }
    else
    {
        // DD channel: B→D X or B→Ds X
        // Parse the two D mesons from the remaining string
        std::string meson1, meson2;

        // Try to match known D meson patterns (order matters - check longer patterns first)
        if (remaining.rfind("dpsdms", 0) == 0)
        {
            meson1 = "dps";
            meson2 = "dms";
        }
        else if (remaining.rfind("dpsd0b", 0) == 0)
        {
            meson1 = "dps";
            meson2 = "d0b";
        }
        else if (remaining.rfind("dpsdm", 0) == 0)
        {
            meson1 = "dps";
            meson2 = "dm";
        }
        else if (remaining.rfind("dpdms", 0) == 0)
        {
            meson1 = "dp";
            meson2 = "dms";
        }
        else if (remaining.rfind("dpd0b", 0) == 0)
        {
            meson1 = "dp";
            meson2 = "d0b";
        }
        else if (remaining.rfind("dpdm", 0) == 0)
        {
            meson1 = "dp";
            meson2 = "dm";
        }
        else
        {
            throw std::runtime_error("Error in parseChannel: Unable to parse DD channel: " + channel);
        }

        return {bMeson, {meson1, meson2}};
    }
}

//----------------------------------------------------------
// Getter for B meson lifetime
double goldenmodesB::getBMesonLifetime(const std::string &bMeson) const
{
    static const std::unordered_map<std::string, double> lifetimes = {
        {"Bp", tau_Bp},
        {"Bd", tau_Bd},
        {"Bs", tau_Bs}};

    auto it = lifetimes.find(bMeson);
    if (it != lifetimes.end())
    {
        return it->second;
    }

    throw std::runtime_error("Error in getBMesonLifetime: Unknown B meson '" + bMeson + "'");
}

//----------------------------------------------------------
// Getter for B meson mass
double goldenmodesB::getBMesonMass(const std::string &bMeson) const
{
    static const std::unordered_map<std::string, double> masses = {
        {"Bp", m_Bp},
        {"Bd", m_Bd},
        {"Bs", m_Bs}};

    auto it = masses.find(bMeson);
    if (it != masses.end())
    {
        return it->second;
    }

    throw std::runtime_error("Error in getBMesonMass: Unknown B meson '" + bMeson + "'");
}

double goldenmodesB::CalculateBR(Parameter amplitude, const std::string &channel) const
{
    // Parse the channel to extract meson components
    auto parsed = parseChannel(channel);
    const std::string &bMeson = parsed.first;
    const std::string &meson1 = parsed.second.first;
    const std::string &meson2 = parsed.second.second;

    // Get the masses of the decaying B meson and final-state mesons
    double m_B = getBMesonMass(bMeson);
    double m1 = mesonMasses.at(meson1);
    double m2 = mesonMasses.at(meson2);

    // Ensure physical kinematics: prevent sqrt of a negative number
    double mass_term1 = (m_B * m_B - (m1 + m2) * (m1 + m2));
    double mass_term2 = (m_B * m_B - (m1 - m2) * (m1 - m2));

    if (mass_term1 < 0 || mass_term2 < 0)
    {
        throw std::runtime_error("Error in CalculateBR: Kinematically forbidden decay for channel " + channel);
    }

    // Compute the magnitude of the final-state momentum
    double p = sqrt(mass_term1 * mass_term2) / (2.0 * m_B);

    // Compute the decay width (Γ)
    double decay_width = (G_F * G_F / (32.0 * M_PI * h_t)) * (p / (m_B * m_B)) * amplitude.Abs2();

    // Compute the branching ratio using the B meson's lifetime
    double lifetime = getBMesonLifetime(bMeson);
    double BR = decay_width * lifetime;

    // Apply symmetry factor for identical particles in the final state
    if (meson1 == meson2)
    {
        BR *= 0.5;
    }

    return BR;
}

// -----------------------------------------------------------------------
// Function to calculate A_CP asymmetry
double goldenmodesB::CalculateAcp(const Parameter &amplitude, const Parameter &conjugate_amplitude) const
{
    // Ensure amplitudes are non-zero to avoid division errors
    double A2 = amplitude.Abs2();
    double Abar2 = conjugate_amplitude.Abs2();

    if (A2 + Abar2 == 0)
    {
        throw std::runtime_error("Error: CalculateAcp - Both amplitudes are zero, division by zero detected.");
    }

    // Compute A_CP
    return (Abar2 - A2) / (A2 + Abar2);
}

// ---------------------------------------------------------------------
// Function to calculate direct CP violation parameter C
double goldenmodesB::CalculateC(const Parameter &amplitude, const Parameter &conjugate_amplitude, const std::string &channel)
{
    // Parse the channel to determine the B meson type
    auto parsed = parseChannel(channel);
    std::string bMeson = parsed.first;

    // Get q/p ratio for Bd or Bs
    TComplex q_p = (bMeson == "Bd") ? ckm.get_q_p_Bd() : ckm.get_q_p_Bs();

    // Special case for K0s and K0l channels (apply q/p_KS)
    if (channel == "Bdjpsik0s" || channel == "Bdjpsik0l" || channel == "Bsjpsik0s")
    {
        q_p *= ckm.get_q_p_KS(); // Multiply by q/p for K0 mixing
    }

    // // Get CP eigenvalue for the channel (ensure it exists)
    // if (cpEigenvalue.find(channel) == cpEigenvalue.end())
    // {
    //     throw std::runtime_error("Error: CalculateC - CP eigenvalue missing for channel: " + channel);
    // }
    // double eta = cpEigenvalue.at(channel);

    // Compute λ (lambda) = η * (q/p) * (A_conjugate / A)
    if (std::abs(amplitude) == 0)
    {
        std::cerr << "Error: CalculateC - Zero amplitude for " << channel << ", division by zero detected." << std::endl;
        return 0.0; // Return safe value instead of crashing
    }

    TComplex lambda = eta * q_p * (conjugate_amplitude / amplitude);

    // Compute C observable: C = (1 - |λ|^2) / (1 + |λ|^2)
    double mod_lambda_squared = lambda.Abs2();
    return (1.0 - mod_lambda_squared) / (1.0 + mod_lambda_squared);
}

// ---------------------------------------------------------------------
// Function to calculate CP violation parameter S
double goldenmodesB::CalculateS(const Parameter &amplitude, const Parameter &conjugate_amplitude, const std::string &channel)
{
    // Parse the channel to determine the B meson type
    auto parsed = parseChannel(channel);
    std::string bMeson = parsed.first;

    // Get q/p ratio for Bd or Bs
    TComplex q_p = (bMeson == "Bd") ? ckm.get_q_p_Bd() : ckm.get_q_p_Bs();

    // Special case for K0s and K0l channels (apply q/p_KS)
    if (channel == "Bdjpsik0s" || channel == "Bdjpsik0l" || channel == "Bsjpsik0s")
    {
        q_p *= ckm.get_q_p_KS(); // Multiply by q/p for K0 mixing
    }

    // // Get CP eigenvalue for the channel (ensure it exists)
    // if (cpEigenvalue.find(channel) == cpEigenvalue.end())
    // {
    //     throw std::runtime_error("Error: CalculateS - CP eigenvalue missing for channel: " + channel);
    // }
    // double eta = cpEigenvalue.at(channel);

    // Compute λ (lambda) = η * (q/p) * (A_conjugate / A)
    if (std::abs(amplitude) == 0)
    {
        std::cerr << "Error: CalculateS - Zero amplitude for " << channel << ", division by zero detected." << std::endl;
        return 0.0; // Return safe value instead of crashing
    }

    TComplex lambda = eta * q_p * (conjugate_amplitude / amplitude);

    // Compute S observable: S = 2 Im(λ) / (1 + |λ|^2)
    double mod_lambda_squared = lambda.Abs2();
    return (2.0 * std::imag(lambda)) / (1.0 + mod_lambda_squared);
}

std::pair<double, double> goldenmodesB::CalculatePhiAndLambda(const Parameter &amplitude, const Parameter &conjugate_amplitude, const std::string &channel)
{
    // Ensure the amplitude is nonzero to avoid division by zero
    if (std::abs(amplitude) == 0)
    {
        throw std::runtime_error("Error: CalculatePhiAndLambda - Zero amplitude detected for channel: " + channel);
    }

    // Get q/p for the Bs meson
    TComplex q_p = ckm.get_q_p_Bs();

    // Compute lambda = (q/p) * (A_cp / A_conj)
    TComplex lambda = q_p * (conjugate_amplitude / amplitude);

    // Compute |lambda|
    double mod_lambda = lambda.Abs();

    // Compute phi_s = -arg(lambda)
    double phi_s = -lambda.Arg();

    return {phi_s, mod_lambda};
}

// std::pair<std::vector<std::string>, std::string> goldenmodesB::extractChannelFromCorrKey(const std::string &corr_key)
// {
//     std::vector<std::string> channels;
//     std::string experiment;

//     if (corr_key.rfind("CS_", 0) == 0)
//     {
//         // Format: "CS_channel_exp"
//         size_t underscore = corr_key.find("_", 3);
//         if (underscore != std::string::npos)
//         {
//             channels.push_back(corr_key.substr(3, underscore - 3)); // Extract channel
//             experiment = corr_key.substr(underscore + 1);           // Extract experiment
//         }
//     }
//     else if (corr_key.rfind("ACP_", 0) == 0)
//     {
//         // Format: "ACP_channel1_channel2_exp"
//         size_t first_underscore = corr_key.find("_", 4);
//         size_t second_underscore = corr_key.find("_", first_underscore + 1);
//         if (first_underscore != std::string::npos && second_underscore != std::string::npos)
//         {
//             channels.push_back(corr_key.substr(4, first_underscore - 4));
//             channels.push_back(corr_key.substr(first_underscore + 1, second_underscore - first_underscore - 1));
//             experiment = corr_key.substr(second_underscore + 1);
//         }
//     }
//     else if (corr_key.rfind("phi_lambda_", 0) == 0)
//     {
//         // Format: "phi_lambda_channel_exp"
//         size_t underscore = corr_key.find("_", 11);
//         if (underscore != std::string::npos)
//         {
//             channels.push_back(corr_key.substr(11, underscore - 11));
//             experiment = corr_key.substr(underscore + 1);
//         }
//     }
//     else if (corr_key.rfind("phi_", 0) == 0)
//     {
//         // Format: "phi_channel_exp"
//         size_t underscore = corr_key.find("_", 4);
//         if (underscore != std::string::npos)
//         {
//             channels.push_back(corr_key.substr(4, underscore - 4));
//             experiment = corr_key.substr(underscore + 1);
//         }
//     }
//     else if (corr_key.rfind("polarization_", 0) == 0)
//     {
//         // Format: "polarization_channel_exp"
//         size_t underscore = corr_key.find("_", 13);
//         if (underscore != std::string::npos)
//         {
//             channels.push_back(corr_key.substr(13, underscore - 13));
//             experiment = corr_key.substr(underscore + 1);
//         }
//     }
//     else
//     {
//         std::cerr << "Warning!! Unknown key format: " << corr_key << std::endl;
//     }

//     return {channels, experiment};
// }

std::map<std::string, double> goldenmodesB::getPolarizationParams(
    const std::string &channel,
    const std::map<std::string, std::pair<Parameter, Parameter>> &amplitude_map)
{
    std::map<std::string, double> polarization_pars;

    try
    {
        // Get parameters for phases
        TComplex B0 = getPar("B_" + channel + "_0");
        TComplex Bperp = getPar("B_" + channel + "_perp");
        TComplex Bparal = getPar("B_" + channel + "_paral");

        double delta_0 = std::arg(B0);
        double deltaparal = std::arg(Bparal);
        double deltaperp = std::arg(Bperp);
        double delta_paral = deltaparal - delta_0;
        double delta_perp = deltaperp - delta_0;

        // Get amplitude information
        auto amp_it = amplitude_map.find(channel);
        if (amp_it == amplitude_map.end())
        {
            throw std::runtime_error("Amplitude map does not contain channel: " + channel);
        }

        TComplex amp = amp_it->second.first;
        TComplex conj_amp = amp_it->second.second;
        TComplex avg_amp = (amp + conj_amp) * 0.5;
        double norm_amp = std::norm(avg_amp);

        if (norm_amp == 0)
        {
            throw std::runtime_error("Normalization factor is zero for channel: " + channel);
        }

        // Get decay amplitudes
        TComplex amp_0 = get_amplitude(channel + "_0");
        TComplex amp_paral = get_amplitude(channel + "_paral");
        TComplex amp_perp = get_amplitude(channel + "_perp");
        TComplex conj_amp_0 = get_conjugate_amplitude(channel + "_0");
        TComplex conj_amp_paral = get_conjugate_amplitude(channel + "_paral");
        TComplex conj_amp_perp = get_conjugate_amplitude(channel + "_perp");

        TComplex avg_A0 = (amp_0 + conj_amp_0) * 0.5;
        TComplex avg_Aparal = (amp_paral + conj_amp_paral) * 0.5;
        TComplex avg_Aperp = (amp_perp + conj_amp_perp) * 0.5;

        double norm_A0 = std::norm(avg_A0);
        double norm_Aparal = std::norm(avg_Aparal);
        double norm_Aperp = std::norm(avg_Aperp);

        // Calculate polarization fractions
        double f_0 = norm_A0 / norm_amp;
        double f_perp = norm_Aperp / norm_amp;
        double f_paral = norm_Aparal / norm_amp;

        // Store the results in the map
        polarization_pars["f_0_" + channel] = f_0;
        polarization_pars["f_perp_" + channel] = f_perp;
        polarization_pars["f_paral_" + channel] = f_paral;
        polarization_pars["delta_paral_" + channel] = delta_paral;
        polarization_pars["delta_perp_" + channel] = delta_perp;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error in getPolarizationParams for channel " << channel << ": " << e.what() << std::endl;
    }

    return polarization_pars;
}

double goldenmodesB::Calculate_UncorrelatedObservables(const std::map<std::string, std::pair<Parameter, Parameter>> &amplitude_map)
{
    double ll_uncorr = 0.0; // Initialize log-likelihood contribution

    for (const auto &channel : channels)
    {
        try
        {
            if (! std::ranges::contains(vectorMesonChannels, channel)) // not a vector meson channel
            {
            auto it = amplitude_map.find(channel);
            if (it == amplitude_map.end())
            {
                std::cerr << "Warning: Amplitude not found for " << channel << std::endl;
                continue;
            }
            const auto &amp_pair = it->second;

            // **Compute ACP, C, and S only if they exist in meas**
            if (meas.find("ACP" + channel) != meas.end())
            {
                double acp = CalculateAcp(amp_pair.first, amp_pair.second);
                obs["ACP_" + channel] = acp;

                double observed = meas.at("ACP" + channel).getMean();
                double uncertainty = meas.at("ACP" + channel).getSigma();
                double diff = acp - observed;
                ll_uncorr += -0.5 * (diff * diff / (uncertainty * uncertainty));
            }

            if (meas.find("C" + channel) != meas.end())
            {
                double c = CalculateC(amp_pair.first, amp_pair.second, channel);
                obs["C_" + channel] = c;

                double observed = meas.at("C" + channel).getMean();
                double uncertainty = meas.at("C" + channel).getSigma();
                double diff = c - observed;
                ll_uncorr += -0.5 * (diff * diff / (uncertainty * uncertainty));
            }

            if (meas.find("S" + channel) != meas.end())
            {
                double s = CalculateS(amp_pair.first, amp_pair.second, channel);
                obs["S_" + channel] = s;

                double observed = meas.at("S" + channel).getMean();
                double uncertainty = meas.at("S" + channel).getSigma();
                double diff = s - observed;
                ll_uncorr += -0.5 * (diff * diff / (uncertainty * uncertainty));
            }

            // **Compute and store branching ratios if they exist in meas**

            std::string br_obsKey;
            std::string br_measKey;
            std::string channel_for_br = channel;
            if (channel == "Bdjpsik0s" || channel == "Bdjpsik0l")
            {
                br_obsKey = "BR_Bdjpsik0";
                br_measKey = "BRBdjpsik0";
                channel_for_br = "Bdjpsik0";             // Use the non-rotated state for BR
                amp_pair = amplitude_map.at("Bdjpsik0"); // Use the non-rotated stated for BR
            }
            else
            {
                br_obsKey = "BR_" + channel;
                br_measKey = "BR" + channel;
            }
            if (meas.find(br_measKey) != meas.end())
            {
                double br_predicted = CalculateBR(amp_pair.first, channel_for_br);
                obs[br_obsKey] = br_predicted;

                double observed = meas.at(br_measKey).getMean();
                double uncertainty = meas.at(br_measKey).getSigma();
                double diff = br_predicted - observed;
                ll_uncorr += -0.5 * (diff * diff / (uncertainty * uncertainty));
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Error in Calculate_UncorrelatedObservables for " << channel << ": " << e.what() << std::endl;
        }
    }

    // **Handle BR Ratios (`R_`)**
    std::vector<std::pair<std::string, std::pair<std::string, std::string>>> br_ratios = {
        {"R_Bpjpsipp_Bpjpsikp", {"Bpjpsipp", "Bpjpsikp"}},
        {"R_Bdjpsiom_Bdjpsirh", {"Bdjpsiom", "Bdjpsirh"}},
        {"R_Bdjpsikst_Bdjpsik0", {"Bdjpsikst", "Bdjpsik0"}},
        {"R_Bdjpsieta_Bsjpsieta", {"Bdjpsieta", "Bsjpsieta"}},
        {"R_Bdjpsietap_Bsjpsietap", {"Bdjpsietap", "Bsjpsietap"}},
        {"R_Bdjpsietap_Bdjpsieta", {"Bdjpsietap", "Bdjpsieta"}},
        {"R_Bsjpsieta_Bsjpsphi", {"Bsjpsieta", "Bsjpsphi"}},
        {"R_Bsjpsietap_Bsjpsieta", {"Bsjpsietap", "Bsjpsieta"}},
        {"R_Bsjpsietap_Bsjpsieta", {"Bsjpsietap", "Bsjpsieta"}},
        {"R_Bsjpsietap_Bsjpsiphi", {"Bsjpsietap", "Bsjpsiphi"}}};

    for (const auto &[ratioKey, channels] : br_ratios)
    {
        if (obs.find("BR_" + channels.first) != obs.end() && obs.find("BR_" + channels.second) != obs.end())
        {
            double BR1 = obs["BR_" + channels.first];
            double BR2 = obs["BR_" + channels.second];

            double R_predicted = BR1 / BR2;
            obs[ratioKey] = R_predicted;

            if (meas.find(ratioKey) != meas.end())
            {
                double observed = meas.at(ratioKey).getMean();
                double uncertainty = meas.at(ratioKey).getSigma();
                double diff = R_predicted - observed;
                ll_uncorr += -0.5 * (diff * diff / (uncertainty * uncertainty));
            }
        }
        else
        {
            std::cerr << "Warning: BR missing for " << ratioKey << std::endl;
        }
    }

    // **Handle Delta A (`deltaA_`)**
    std::vector<std::pair<std::string, std::pair<std::string, std::string>>> deltaA_ratios = {
        {"deltaA_Bpjpsipp_Bpjpsikp", {"ACP_Bpjpsipp", "ACP_Bpjpsikp"}}};

    for (const auto &[deltaAKey, channels] : deltaA_ratios)
    {
        if (obs.find(channels.first) != obs.end() && obs.find(channels.second) != obs.end())
        {
            double ACP1 = obs[channels.first];
            double ACP2 = obs[channels.second];

            double deltaA_predicted = ACP1 - ACP2;
            obs[deltaAKey] = deltaA_predicted;

            if (meas.find(deltaAKey) != meas.end())
            {
                double observed = meas.at(deltaAKey).getMean();
                double uncertainty = meas.at(deltaAKey).getSigma();
                double diff = deltaA_predicted - observed;
                ll_uncorr += -0.5 * (diff * diff / (uncertainty * uncertainty));
            }
        }
        else
        {
            std::cerr << "Warning: ACP missing for " << deltaAKey << std::endl;
        }
    }

    return ll_uncorr;
}

//----------------------------------------------------------------------------------

double goldenmodesB::Calculate_CorrelatedObservables(const std::map<std::string, std::pair<Parameter, Parameter>> &amplitude_map)
{
    double ll_corr = 0.0; // Initialize the log-likelihood contribution
    TVectorD corr(2);     // For correlated observables (e.g., C and S)
    TVectorD corr4(4);
    TVectorD corr5(5);
    TVectorD corr6(6);

    // **Correlated observables for Bdjpsik0s**
    {
        const auto &amp_pair_Bdjpsik0s = amplitude_map.at("Bdjpsik0s");
        if (obs.find("C_Bdjpsik0s") != obs.end() && obs.find("S_Bdjpsik0s") != obs.end())
        {
            corr(0) = obs["C_Bdjpsik0s"];
            corr(1) = obs["S_Bdjpsik0s"];

            ll_corr += corrmeas.at("CS_Bdjpsik0s_LHCb2023").logweight(corr);
            ll_corr += corrmeas.at("CS_Bdjpsik0s_BelleII2024").logweight(corr);
        }
        else
        {
            std::cerr << "Error: C and S values for Bdjpsik0s not found in obs map!" << std::endl;
        }
    }

    // **Correlated observables for Bdjpsip0**
    {
        const auto &amp_pair_Bdjpsip0 = amplitude_map.at("Bdjpsip0");
        if (obs.find("C_Bdjpsip0") != obs.end() && obs.find("S_Bdjpsip0") != obs.end())
        {
            corr(0) = obs["C_Bdjpsip0"];
            corr(1) = obs["S_Bdjpsip0"];
            ll_corr += corrmeas.at("CS_Bdjpsip0_BaBar2008").logweight(corr);
        }
        else
        {
            std::cerr << "Error: C and S values for Bdjpsip0 not found in obs map!" << std::endl;
        }
    }

    // **Correlated observables for Bdjpsirh**
    {
        try
        {
            const auto &amp_pair_Bdjpsirh = amplitude_map.at("Bdjpsirh");
            double c_predicted = CalculateC(amp_pair_Bdjpsirh.first, amp_pair_Bdjpsirh.second, "Bdjpsirh");
            double s_predicted = CalculateS(amp_pair_Bdjpsirh.first, amp_pair_Bdjpsirh.second, "Bdjpsirh");
            obs["C_Bdjpsirh"] = c_predicted;
            obs["S_Bdjpsirh"] = s_predicted;
            corr(0) = c_predicted;
            corr(1) = s_predicted;
            ll_corr += corrmeas.at("CS_Bdjpsirh_LHCb2014").logweight(corr);
        }
        catch (const std::out_of_range &e)
        { // Catch exception for missing key
            std::cerr << "Error: Bdjpsirh not found in amplitude_map! " << e.what() << std::endl;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Unexpected error in correlated observables for Bdjpsirh: " << e.what() << std::endl;
        }
    }

    // **Correlated observables for Bsjpsiph**
    {
        auto pol_params = getPolarizationParams("Bsjpsiph", amplitude_map);

        // Extract the polarization fractions and phases:
        double f_0 = pol_params.at("f_0_Bsjpsiph");
        double f_paral = pol_params.at("f_paral_Bsjpsiph");
        double f_perp = pol_params.at("f_perp_Bsjpsiph");
        double delta_paral = pol_params.at("delta_paral_Bsjpsiph");
        double delta_perp = pol_params.at("delta_perp_Bsjpsiph");

        obs["f_0_Bsjpsiph"] = f_0;
        obs["f_paral_Bsjpsiph"] = f_paral;
        obs["f_perp_Bsjpsiph"] = f_perp;
        obs["delta_paral_Bsjpsiph"] = delta_paral;
        obs["delta_perp_Bsjpsiph"] = delta_perp;

        // --- Retrieve additional parameters (e.g. weak phase and |lambda|) ---
        double phi_Bsjpsiph = obs["phi_Bsjpsiph"];
        double lambda_Bsjpsiph = obs["lambda_Bsjpsiph"];

        // Loop over the correlated measurements
        for (auto &[key, corrObs] : corrmeas)
        {
            if (key.find("Bsjpsiph") != std::string::npos)
            {
                // Compute phi and lambda
                const auto &amp_pair_Bsjpsiph = amplitude_map.at("Bsjpsiph");
                // auto phi_lambda_result = CalculatePhiAndLambda(amp_pair_Bsjpsiph.first, amp_pair_Bsjpsiph.second, "Bsjpsiph");
                double phi_Bsjpsiph = obs["phi_Bsjpsiph"];
                double lambda_Bsjpsiph = obs["lambda_Bsjpsiph"];
                // Handle each correlated measurement key
                if (key == "phi_lambda_Bsjpsiph_LHCb2023")
                {
                    TVectorD corr6(6);          // 6 correlated observables
                    corr6(0) = phi_Bsjpsiph;    // Weak phase
                    corr6(1) = lambda_Bsjpsiph; // |lambda|
                    corr6(2) = f_perp;          // Polarization fraction
                    corr6(3) = f_0;             // Polarization fraction
                    corr6(4) = delta_perp;      // Strong phase
                    corr6(5) = delta_paral;     // Strong phase
                    ll_corr += corrObs.logweight(corr6);
                }
                else if (key == "phi_Bsjpsiph_ATLAS2020B")
                {
                    TVectorD corr5(5); // 5 correlated observables
                    corr5(0) = phi_Bsjpsiph;
                    corr5(1) = f_paral;
                    corr5(2) = f_0;
                    corr5(3) = delta_perp;
                    corr5(4) = delta_paral;
                    ll_corr += corrObs.logweight(corr5);
                }
                else if (key == "phi_lambda_Bsjpsiph_LHCb2021")
                {
                    TVectorD corr6(6); // 6 correlated observables
                    corr6(0) = phi_Bsjpsiph;
                    corr6(1) = lambda_Bsjpsiph;
                    corr6(2) = f_perp;
                    corr6(3) = f_0;
                    corr6(4) = delta_perp;
                    corr6(5) = delta_paral;
                    ll_corr += corrObs.logweight(corr6);
                }
                else if (key == "phi_Bsjpsiph_CMS2020")
                {
                    TVectorD corr5(5); // 5 correlated observables
                    corr5(0) = phi_Bsjpsiph;
                    corr5(1) = f_perp;
                    corr5(2) = f_0;
                    corr5(3) = delta_perp;
                    corr5(4) = delta_paral;
                    ll_corr += corrObs.logweight(corr5);
                }
                else if (key == "phi_lambda_Bsjpsiph_LHCb2019" || key == "phi_lambda_Bsjpsiph_LHCb2017")
                {
                    TVectorD corr2(2); // Only phi and lambda
                    corr2(0) = phi_Bsjpsiph;
                    corr2(1) = lambda_Bsjpsiph;
                    ll_corr += corrObs.logweight(corr2);
                }
                else if (key == "phi_Bsjpsiph_ATLAS2016")
                {
                    TVectorD corr5(5); // 5 correlated observables
                    corr5(0) = phi_Bsjpsiph;
                    corr5(1) = f_paral;
                    corr5(2) = f_0;
                    corr5(3) = delta_perp;
                    corr5(4) = delta_paral;
                    ll_corr += corrObs.logweight(corr5);
                }
                else if (key == "phi_Bsjpsiph_ATLAS2014")
                {
                    TVectorD corr5(5); // 5 correlated observables
                    corr5(0) = phi_Bsjpsiph;
                    corr5(1) = f_paral;
                    corr5(2) = f_0;
                    corr5(3) = delta_perp;
                    corr5(4) = delta_paral;
                    ll_corr += corrObs.logweight(corr5);
                }
            }
        }
    }

    // correlated polarization measurments for C_Bdjpsikst
    {
        auto pol_params = getPolarizationParams("Bdjpsikst", amplitude_map);

        // Extract the polarization fractions and phases:
        double f_0 = pol_params.at("f_0_Bdjpsikst");
        double f_paral = pol_params.at("f_paral_Bdjpsikst");
        double f_perp = pol_params.at("f_perp_Bdjpsikst");
        double delta_paral = pol_params.at("delta_paral_Bdjpsikst");
        double delta_perp = pol_params.at("delta_perp_Bdjpsikst");

        obs["f_0_Bdjpsikst"] = f_0;
        obs["f_paral_Bdjpsikst"] = f_paral;
        obs["f_perp_Bdjpsikst"] = f_perp;
        obs["delta_paral_Bdjpsikst"] = delta_paral;
        obs["delta_perp_Bdjpsikst"] = delta_perp;

        TVectorD corr4(4);      // 4 correlated observables
        corr4(0) = f_paral;     // Polarization fraction
        corr4(1) = f_perp;      // Polarization fraction
        corr4(2) = delta_paral; // Strong phase
        corr4(3) = delta_perp;  // Strong phase
        ll_corr += corrmeas.at("polarization_Bdjpsikst_LHCb2013").logweight(corr4);
    }

    {
        auto pol_params = getPolarizationParams("Bsjpsikst", amplitude_map);

        // Extract the polarization fractions and phases:
        double f_0 = pol_params.at("f_0_Bsjpsikst");
        double f_paral = pol_params.at("f_paral_Bsjpsikst");
        double delta_paral = pol_params.at("delta_paral_Bsjpsikst");
        double delta_perp = pol_params.at("delta_perp_Bsjpsikst");

        obs["f_0_Bsjpsikst"] = f_0;
        obs["f_paral_Bsjpsikst"] = f_paral;
        obs["delta_paral_Bsjpsikst"] = delta_paral;
        obs["delta_perp_Bsjpsikst"] = delta_perp;
        // --- Retrieve additional parameters  ---
        complex<double> amp_0 = get_amplitude("Bsjpsikst_0");
        complex<double> amp_paral = get_amplitude("Bsjpsikst_paral");
        complex<double> amp_perp = get_amplitude("Bsjpsikst_perp");
        complex<double> conj_amp_0 = get_conjugate_amplitude("Bsjpsikst_0");
        complex<double> conj_amp_paral = get_conjugate_amplitude("Bsjpsikst_paral");
        complex<double> conj_amp_perp = get_conjugate_amplitude("Bsjpsikst_perp");

        obs["A0_CP_Bsjpsikst"] = CalculateAcp(amp_0, conj_amp_0);
        obs["Aparal_CP_Bsjpsikst"] = CalculateAcp(amp_paral, conj_amp_paral);
        obs["Aperp_CP_Bsjpsikst"] = CalculateAcp(amp_perp, conj_amp_perp);

        TVectorD corr7(7);
        corr7(0) = obs["A0_CP_Bsjpsikst"];
        corr7(1) = obs["Aparal_CP_Bsjpsikst"];
        corr7(2) = obs["Aperp_CP_Bsjpsikst"];
        corr7(3) = f_0;
        corr7(4) = f_paral;
        corr7(5) = delta_paral;
        corr7(6) = delta_perp;
        ll_corr += corrmeas.at("polarization_Bsjpsikst_LHCb2015").logweight(corr7);
    }
    return ll_corr;
}

//------------------------------------------------------------

double goldenmodesB::LogLikelihood(const std::vector<double> &parameters)
{
    static int iteration_counter = 0;
    ++iteration_counter;

    obs.clear();     // Clear obs map for each iteration
    double ll = 0.0; // Log-likelihood accumulator

    // Check that the parameter vector has the expected size
    int expectedSize = 0;
    for (const auto &channel : channelNamesSU3)
    {
        auto it = channelParameters.find(channel);
        if (it != channelParameters.end())
        {
            expectedSize += it->second.size();
        }
    }
    // Last 5 parameters are for CKM + eta-eta' mixing angle
    if (parameters.size() != expectedSize + 5)
    {
        std::cerr << "Error: parameters.size() = " << parameters.size()
                  << ", but expected size = " << expectedSize + 5 << std::endl;
        return log(0.); // Return log(0) for invalid size
    }
    if (parameters.empty())
    {
        std::cerr << "Error: Empty parameters vector!" << std::endl;
        return log(0.);
    }

    // Unpack CKM parameters and compute CKM elements
    std::vector<double> ckmParams(parameters.end() - 5, parameters.end() - 1);
    ckm.computeCKM(ckmParams[0], ckmParams[1], ckmParams[2], ckmParams[3], true);

    // eta-eta' mixing angle
    double theta_P = parameters.back();
    SetParameterValue("theta_P", theta_P);

    // Unpack the parameters for each channel (from channelNamesSU3)
    int index = 0;
    for (const auto &channel : channelNamesSU3)
    {
        auto it = channelParameters.find(channel);
        if (it != channelParameters.end())
        {
            // Loop over parameter pairs (real and imaginary)
            for (size_t j = 0; j < it->second.size(); j += 2)
            {
                if (index + 1 >= expectedSize)
                { // Only loop over channel parameters
                    break;
                }

                double realPart = parameters[index];
                double imagPart = parameters[index + 1];
                const std::string &realParamName = it->second[j];
                const std::string &imagParamName = it->second[j + 1];

                SetParameterValue(realParamName, realPart);
                SetParameterValue(imagParamName, imagPart);

                index += 2;
            }
        }
        else
        {
            std::cerr << "Channel " << channel << " not found in channelParameters." << std::endl;
            return log(0.); // Return log(0) for invalid channel
        }
    }
    // Print parameters every 1000 iterations
    /*if (iteration_counter % 1000 == 0) {
        std::cout << "Iteration " << iteration_counter << " - Parameter Values:" << std::endl;
        for (const auto& param : parameterValues) {
            std::cout << param.first << " = " << param.second << std::endl;
        }
    }*/
    // Build the amplitude map using physical amplitudes

    amplitude_map.clear();

    for (const std::string &channel : channelNamesSU3)
    {
        if (std::ranges::contains(vectorMesonChannels, channel))
        {
            amplitude_map[channel + "_0"] = {get_amplitude(channel + "_0"), get_conjugate_amplitude(channel + "_0")};
            amplitude_map[channel + "_perp"] = {get_amplitude(channel + "_perp"), get_conjugate_amplitude(channel + "_perp")};
            amplitude_map[channel + "_paral"] = {get_amplitude(channel + "_paral"), get_conjugate_amplitude(channel + "_paral")};
        }
        else
        {
            // Directly store original amplitudes
            amplitude_map[channel] = {get_amplitude(channel), get_conjugate_amplitude(channel)};

            // For CP observables, add K0s and K0l cases with the same amplitudes as Bdjpsik0
            if (channel == "Bdjpsik0")
            {
                amplitude_map["Bdjpsik0s"] = amplitude_map["Bdjpsik0"];
                amplitude_map["Bdjpsik0l"] = amplitude_map["Bdjpsik0"];
            }
            else if (channel == "Bsjpsik0")
            {
                amplitude_map["Bsjpsik0s"] = amplitude_map["Bsjpsik0"];
            }
        }

        for (const std::string &channel : channelNames)
        {
            // compute amplitudes for eta and etaprime final states
            if (channel.ends_with("eta"))
            {
                amplitude_map[channel] = {cos(theta_P) * amplitude_map[channel + "8"].first - sin(theta_P) * amplitude_map[channel + "1"].first,
                                          cos(theta_P) * amplitude_map[channel + "8"].second - sin(theta_P) * amplitude_map[channel + "1"].second};
                amplitude_map[channel + "p"] = {sin(theta_P) * amplitude_map[channel + "8"].first + cos(theta_P) * amplitude_map[channel + "1"].first,
                                                sin(theta_P) * amplitude_map[channel + "8"].second + cos(theta_P) * amplitude_map[channel + "1"].second};
            }
        }

        // Check NaN/Inf values in amplitudes
        for (const auto &[chan, amp_pair] : amplitude_map)
        {
            if (std::isnan(abs(amp_pair.first)) || std::isinf(abs(amp_pair.first)) ||
                std::isnan(abs(amp_pair.second)) || std::isinf(abs(amp_pair.second)))
            {
                std::cerr << "Invalid amplitude (NaN or Inf) for channel: " << chan << std::endl;
                return log(0.);
            }
        }

        // Add contributions from uncorrelated and correlated observables
        ll += Calculate_UncorrelatedObservables(amplitude_map);
        ll += Calculate_CorrelatedObservables(amplitude_map);

        return ll;
    }

    //---------------------------------------------------------

    void goldenmodesB::MCMCUserIterationInterface()
    {
        // Loop over all MCMC chains
        const unsigned int log_interval = 1000; // Log every 1000 iterations
        std::vector<double> pars;

        for (unsigned int i = 0; i < fMCMCNChains; ++i)
        {
            pars = fMCMCStates.at(i).parameters;
            try
            {
                LogLikelihood(pars);
            }
            catch (const std::exception &e)
            {
                std::cerr << "Error in LogLikelihood: " << e.what() << std::endl;
                return;
            }

            // Directly fill histograms
            histos.fillh1d();
            histos.fillh2d();
        }
    }

    void goldenmodesB::SaveHistograms(const std::string &filename)
    {
        TFile file(filename.c_str(), "RECREATE");

        // Write all histograms to the file
        for (const auto &histPair : histos.h1d)
        {
            histPair.second->Write();
        }

        // Write all 2D histograms
        for (const auto &histPair : histos.h2d)
        {
            histPair.second->Write();
        }

        file.Close();
        std::cout << "Histograms saved to " << filename << std::endl;
    }

    // ---------------------------------------------------------

    void goldenmodesB::PrintObservablePulls(const std::string &filename)
    {
        std::ofstream outfile(filename);

        if (!outfile.is_open())
        {
            std::cerr << "Error: Unable to open file " << filename << std::endl;
            return;
        }

        // Header for the output file
        outfile << "Observable\t\t Measurement\t\t Mean\t\t Sigma\t\t Pull\n";

        // Define mappings for observable names and measurement keys
        std::map<std::string, std::string> measMap = {
            {"BR_Bdjpsik0", "BRBdjpsik0"}, {"BR_Bdjpsip0", "BRBdjpsip0"}, {"BR_Bdjpsiom", "BRBdjpsiom"}, {"BR_Bpjpsikp", "BRBpjpsikp"}, {"BR_Bpjpsipp", "BRBpjpsipp"}, {"BR_Bsjpsiph", "BRBsjpsiph"}, {"BR_Bsjpsik0s", "BRBsjpsik0s"}, {"BR_Bdjpsirh", "BRBdjpsirh"}, {"BR_Bdjpsikst", "BRBdjpsikst"}, {"BR_Bsjpsikst", "BRBsjpsikst"}, {"C_Bdjpsik0s", "CBdjpsik0s"}, {"C_Bdjpsik0l", "CBdjpsik0l"}, {"S_Bdjpsik0s", "SBdjpsik0s"}, {"S_Bdjpsik0l", "SBdjpsik0l"}, {"C_Bdjpsip0", "CBdjpsip0"}, {"S_Bdjpsip0", "SBdjpsip0"}, {"ACP_Bpjpsikp", "ACPBpjpsikp"}, {"ACP_Bpjpsipp", "ACPBpjpsipp"}, {"deltaA_Bpjpsipp_Bpjpsikp", "deltaA_Bpjpsipp_Bpjpsikp"}, {"R_Bpjpsipp_Bpjpsikp", "R_Bpjpsipp_Bpjpsikp"}, {"lambda_Bsjpsiph", "lambda_Bsjpsiph"}, {"phi_Bsjpsiph", "phi_Bsjpsiph"}, {"f_perp_Bsjpsiph", "f_perp_Bsjpsiph"}, {"f_paral_Bsjpsiph", "f_paral_Bsjpsiph"}, {"f_0_Bsjpsiph", "f_0_Bsjpsiph"}, {"delta_perp_Bsjpsiph", "delta_perp_Bsjpsiph"}, {"delta_paral_Bsjpsiph", "delta_paral_Bsjpsiph"}, {"R_Bdjpsiom_Bdjpsirh", "R_Bdjpsiom_Bdjpsirh"}, {"R_Bdjpsikst_Bdjpsik0", "R_Bdjpsikst_Bdjpsik0"}, {"C_Bdjpsikst", "CBdjpsikst"}, {"S_Bdjpsikst", "SBdjpsikst"}, {"f_0_Bdjpsikst", "f_0_Bdjpsikst"}, {"f_paral_Bdjpsikst", "f_paral_Bdjpsikst"}, {"f_perp_Bdjpsikst", "f_perp_Bdjpsikst"}, {"delta_paral_Bdjpsikst", "delta_paral_Bdjpsikst"}, {"delta_perp_Bdjpsikst", "delta_perp_Bdjpsikst"}, {"f_0_Bsjpsikst", "f_0_Bsjpsikst"}, {"f_paral_Bsjpsikst", "f_paral_Bsjpsikst"}, {"delta_paral_Bsjpsikst", "delta_paral_Bsjpsikst"}, {"delta_perp_Bsjpsikst", "delta_perp_Bsjpsikst"}, {"A0_CP_Bsjpsikst", "A0_CP_Bsjpsikst"}, {"Aparal_CP_Bsjpsikst", "Aparal_CP_Bsjpsikst"}, {"Aperp_CP_Bsjpsikst", "Aperp_CP_Bsjpsikst"}, {"BR_Bddpdm", "BRBddpdm"}, {"C_Bddpdm", "CBddpdm"}, {"S_Bddpdm", "SBddpdm"}, {"BR_Bsdpsdms", "BRBsdpsdms"}, {"C_Bsdpsdms", "CBsdpsdms"}, {"S_Bsdpsdms", "SBsdpsdms"}, {"BR_Bpdpd0b", "BRBpdpd0b"}, {"ACP_Bpdpd0b", "ACPBpdpd0b"}, {"BR_Bpdpsd0b", "BRBpdpsd0b"}, {"ACP_Bpdpsd0b", "ACPBpdpsd0b"}

        };

        for (const auto &histPair : histos.h1d)
        {
            const std::string &obs_name = histPair.first;
            TH1D *hist = histPair.second;

            double obs_mean = hist->GetMean();
            double obs_measurement = 0.0;
            double sigma_measurement = 0.0;
            bool found_measurement = false;

            // Check if the observable exists in `newmeas` first
            if (measMap.find(obs_name) != measMap.end())
            {
                const std::string &key = measMap[obs_name];

                if (newmeas.find(key) != newmeas.end())
                {
                    obs_measurement = newmeas.at(key).getMean();
                    sigma_measurement = newmeas.at(key).getSigma();
                    found_measurement = true;
                }
                else if (meas.find(key) != meas.end())
                {
                    // Use `meas` only if it is not found in `newmeas`
                    obs_measurement = meas.at(key).getMean();
                    sigma_measurement = meas.at(key).getSigma();
                    found_measurement = true;
                }
            }

            if (found_measurement)
            {
                // Calculate the pull value: (predicted - observed) / uncertainty
                double pull = (obs_mean - obs_measurement) / sigma_measurement;

                outfile << obs_name << "\t"
                        << obs_measurement << "\t"
                        << obs_mean << "\t"
                        << sigma_measurement << "\t"
                        << pull << "\n";
            }
        }

        outfile.close();
        std::cout << "Pull values saved to " << filename << std::endl;
    }
