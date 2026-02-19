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
#include <regex>
#include <TRandom3.h>
#include <cmath>
#include <TMatrixD.h>
#include <TDecompChol.h>
#include <iostream>
#include <fstream>
using namespace std;

// ---------------------------------------------------------

CKMParameters ckm;

goldenmodesB::goldenmodesB(double &dsu3_limit_in, double &ewp_limit_in, bool BJPSIP, bool BJPSIV, bool BDDb, bool Debug_in) : BCModel(), histos(obs), Debug(Debug_in)
{
    TH1::SetDefaultBufferSize(1000000);
    dsu3_limit = dsu3_limit_in;
    ewp_limit = ewp_limit_in;
    cout << "constructor for goldenmodes called: inserting the experimental data" << endl;
    vector<dato> data;
    PDGAverage pdgaverage;
    evaluatedevts = 0;

    if (BJPSIP)
    {
        channels.insert(channels.end(), pseudoscalarMesonChannels.begin(), pseudoscalarMesonChannels.end());
        channelNamesSU3.insert(channelNamesSU3.end(), pseudoscalarMesonChannelsSU3.begin(), pseudoscalarMesonChannelsSU3.end());
    }
    if (BJPSIV)
    {
        channels.insert(channels.end(), vectorMesonChannels.begin(), vectorMesonChannels.end());
        channelNamesSU3.insert(channelNamesSU3.end(), vectorMesonChannelsSU3.begin(), vectorMesonChannelsSU3.end());
        vector<string> polarizedChannels;
        for (const auto &chan : vectorMesonChannels)
        {
            polarizedChannels.push_back(chan + "_0");
            polarizedChannels.push_back(chan + "_paral");
            polarizedChannels.push_back(chan + "_perp");
        }
        channels.insert(channels.end(), polarizedChannels.begin(), polarizedChannels.end());
    }
    if (BDDb)
    {
        channels.insert(channels.end(), ddbarChannels.begin(), ddbarChannels.end());
        channelNamesSU3.insert(channelNamesSU3.end(), ddbarChannelsSU3.begin(), ddbarChannelsSU3.end());
    }

    DeclareParameters(); // Ensure parameters are defined

    // Add delta phi parameters for Bd and Bs mixing to allow for NP contributions and/or to decorrelate from CKM fit
    AddParameter("myphid", -M_PI / 2., M_PI / 2.); // in rad
    AddParameter("myphis", -M_PI / 2., M_PI / 2.); // in rad

    // Add CKM parameters directly in the constructor
    AddParameter("CKM_Vud", 0.97432 - 5 * 0.00015, 0.97432 + 5 * 0.00015);                       // Vud parameter - 5 sigma range
    AddParameter("CKM_Vcb", 0.04118 - 5 * 0.00076, 0.04118 + 5 * 0.00076);                       // Vcb parameter - 5 sigma range
    AddParameter("CKM_Vub", 0.00382 - 5 * 0.00034, 0.00382 + 5 * 0.00034);                       // Vub parameter - 5 sigma range
    AddParameter("CKM_gamma", (65.7 - 5 * 2.5) / 180.0 * M_PI, (65.7 + 5 * 2.5) / 180.0 * M_PI); // gamma parameter (in rad) - 5 sigma range

    // Add mixing angle between eta1 and eta8
    AddParameter("theta_P", -30. / 180.0 * M_PI, 0.); // in rad

    SetPriorConstantAll();

    // Use lattice QCD result for theta_P
    GetParameter("theta_P").SetPrior(make_shared<BCGaussianPrior>(-14.1 / 180.0 * M_PI, 2.8 / 180.0 * M_PI)); // in rad

    // measurements
    vector<dato> CorrData;
    TMatrixDSym CorrStat(2);
    TMatrixDSym CorrSyst(2);
    vector<string> names;

    if (BJPSIP)
    {

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
        corrmeas_channels.insert(pair<string, vector<string>>("CS_Bdjpsik0s_LHCb2023", names));
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
        corrmeas_channels.insert(pair<string, vector<string>>("CS_Bdjpsik0s_BelleII2024", names));
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
        data.push_back(dato(0.59, 0.19, 0.03)); ////Belle:2018nxw
        data.push_back(dato(0.88, 0.17, 0.03)); // Belle-II:2024hqw

        pdgaverage.setData(data);
        pdgaverage.setName("SBdjpsip0");
        pdgaverage.CalculateAverage();

        meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
        data.clear();

        // Bdjpsip0 : C and S observables from BaBar:2008kfx
        CorrData.push_back(dato(-0.20, 0.19, 0.03)); // C observable
        names.push_back("CBdjpsip0");
        CorrData.push_back(dato(1.23, 0.21, 0.04)); // S observable
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

        corrmeas_channels.insert(pair<string, vector<string>>("CS_Bdjpsip0_BaBar2008", names));
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
        meas.insert(pair<string, dato>("R_Bdjpsieta_Bsjpsieta", dato(2.16e-2, 0.16e-2, 0.05e-2, 0.07e-2))); // LHCb:2025sgp supersedes previous LHCb:2014oms

        ////////////////////////////
        // Bdjpsieta'
        ////////////////////////////

        // Ratios of BRs

        // BR(Bd->J/psi eta')/BR(Bs->J/psi eta')
        meas.insert(pair<string, dato>("R_Bdjpsietap_Bsjpsietap", dato(1.33e-2, 0.12e-2, 0.05e-2, 0.04e-2))); // LHCb:2025sgp supersedes previous LHCb:2014oms

        // BR(Bd->J/psi eta')/BR(Bd->J/psi eta)
        meas.insert(pair<string, dato>("R_Bdjpsietap_Bdjpsieta", dato(0.48, 0.06, 0.02, 0.01))); // LHCb:2025sgp supersedes previous LHCb:2014oms

        ////////////////////////////
        // Bpjpsikp
        ////////////////////////////

        // BR measurements
        data.push_back(dato(10.32e-4, 0.07e-4, 0.24e-4)); // BELLE:2019xld
        data.push_back(dato(9.4e-4, 0.7e-4, 0.8e-4));     // Belle:2019avj
        data.push_back(dato(8.9e-4, 0.6e-4, 0.5e-4));     // Belle:2017psv
        data.push_back(dato(8.1e-4, 1.3e-4, 0.7e-4));     // BaBar:2005pcw
        data.push_back(dato(10.61e-4, 0.15e-4, 0.48e-4)); // BaBar:2004htr
        data.push_back(dato(10.4e-4, 1.1e-4, 0.1e-4));    // BaBar:2005sdl
        data.push_back(dato(10.2e-4, 0.8e-4, 0.7e-4));    // CLEO:1997ilq
        data.push_back(dato(9.24e-4, 3.04e-4, 0.05e-4));  // CLEO:1991roe
        data.push_back(dato(8.09e-4, 3.50e-4, 0.04e-4));  // ARGUS:1990jet

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
        CorrStat(0, 1) = -0.06; // Correlation
        CorrStat(1, 0) = -0.06; // Symmetric part (same as Corr(0, 1))
        CorrSyst.ResizeTo(2, 2);
        CorrSyst.UnitMatrix(); // No systematic correlation available in the paper
        // Insert correlated data into corrmeas
        corrmeas.insert(pair<string, CorrelatedGaussianObservables>("CS_Bsjpsik0s_LHCb2015", CorrelatedGaussianObservables(CorrData, CorrStat, CorrSyst)));
        corrmeas_channels.insert(pair<string, vector<string>>("CS_Bsjpsik0s_LHCb2015", names));
        CorrData.clear();
        names.clear();

        /////////////////////////////
        // Bsjpsieta
        /////////////////////////////

        // BR measurements

        // BRBsjpsieta
        meas.insert(pair<string, dato>("BRBsjpsieta", dato(5.27e-4, 0.50e-4, 0.25e-4, 0.97e-4))); // Chang:2012gnb with symmetrized errors

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
    }
    if (BJPSIP && BJPSIV)
    {

        // R_Bsjpsietap_Bsjpsiphi
        meas.insert(pair<string, dato>("R_Bsjpsietap_Bsjpsiphi", dato(0.370, 0.013, 0.018, 0.011)));        // LHCb:2025sgp
        meas.insert(pair<string, dato>("R_Bsjpsieta_Bsjpsiphi", dato(4.50e-1, 0.14e-1, 0.16e-1, 0.13e-1))); // LHCb:2025sgp
    }

    if (BJPSIV)
    {

        /////////////////////////////
        // Bsjpsiphi
        /////////////////////////////

        // BRBsjpsiphi measurements
        data.push_back(dato(1.018e-3, 0.032e-3, 0.037e-3)); // LHCb:2021qbv
        data.push_back(dato(1.25e-3, 0.07e-3, 0.23e-3));    // Belle:2013sdi
        data.push_back(dato(1.5e-3, 0.5e-3, 0.1e-3));       // CDF:1996ivk

        pdgaverage.setData(data);
        pdgaverage.setName("BRBsjpsiphi");
        pdgaverage.CalculateAverage();

        meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
        data.clear();

        // Correlated data sets for polarization fractions and CP asymmetries

        // LHCb:2023sim (Run 1+2 combined)

        CorrData.push_back(dato(-0.044, 0.020, 0., 0., 0., true)); // phi_s
        names.push_back("phis_Bsjpsiphi");
        CorrData.push_back(dato(0.990, 0.010)); // lambda
        names.push_back("lambda_Bsjpsiphi");
        CorrData.push_back(dato(0.2471, 0.0032)); // f_perp
        names.push_back("f_perp_Bsjpsiphi");
        CorrData.push_back(dato(0.5175, 0.0035)); // f_0
        names.push_back("f_0_Bsjpsiphi");
        CorrData.push_back(dato(2.924, 0.076, 0., 0., 0., true)); // delta_perp - delta_0
        names.push_back("delta_perp_Bsjpsiphi");
        CorrData.push_back(dato(3.150, 0.062, 0., 0., 0., true)); // delta_paral - delta_0
        names.push_back("delta_paral_Bsjpsiphi");

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
        corrmeas.insert(pair<string, CorrelatedGaussianObservables>("Bsjpsiphi_LHCb2023sim", CorrelatedGaussianObservables(CorrData, CorrPhis)));
        corrmeas_channels.insert(pair<string, vector<string>>("Bsjpsiphi_LHCb2023sim", names));
        CorrData.clear();
        names.clear();

        // CMS:2024znt (Run 1+2 combined)

        CorrData.push_back(dato(-0.074, 0.023, 0., 0., 0., true)); // phi_s
        names.push_back("phis_Bsjpsiphi");
        CorrData.push_back(dato(1.011, 0.019)); // lambda
        names.push_back("lambda_Bsjpsiphi");
        CorrData.push_back(dato(0.5273, 0.0044)); // f_0
        names.push_back("f_0_Bsjpsiphi");
        CorrData.push_back(dato(0.2417, 0.0036)); // f_perp
        names.push_back("f_perp_Bsjpsiphi");
        CorrData.push_back(dato(3.152, 0.077, 0. , 0., 0., true)); // delta_paral - delta_0
        names.push_back("delta_paral_Bsjpsiphi");
        CorrData.push_back(dato(2.940, 0.098, 0. , 0., 0., true)); // delta_perp - delta_0
        names.push_back("delta_perp_Bsjpsiphi");

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
        corrmeas.insert(pair<string, CorrelatedGaussianObservables>("Bsjpsiphi_CMS2024znt", CorrelatedGaussianObservables(CorrData, CorrPhis)));
        corrmeas_channels.insert(pair<string, vector<string>>("Bsjpsiphi_CMS2024znt", names));
        CorrData.clear();
        names.clear();

        // phi and lambda for ATLAS:2020lbz (Run 1+2 combined)
        TMatrixDSym Corr5(5);

        CorrData.push_back(dato(-0.087, 0.036, 0.021, 0., 0., true)); // phi_s
        names.push_back("phis_Bsjpsiphi");
        CorrData.push_back(dato(0.2218, 0.0017, 0.0021)); // f_paral
        names.push_back("f_paral_Bsjpsiphi");
        CorrData.push_back(dato(0.5152, 0.0012, 0.0034)); // f_0
        names.push_back("f_0_Bsjpsiphi");
        CorrData.push_back(dato(3.03, 0.10, 0.05, 0., 0., true)); // delta_perp - delta_0
        names.push_back("delta_perp_Bsjpsiphi");
        CorrData.push_back(dato(2.95, 0.05, 0.09, 0., 0., true)); // delta_par - delta_0
        names.push_back("delta_paral_Bsjpsiphi");
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

        corrmeas.insert(pair<string, CorrelatedGaussianObservables>("phi_Bsjpsiphi_ATLAS2020B", CorrelatedGaussianObservables(CorrData, Corr5)));
        corrmeas_channels.insert(pair<string, vector<string>>("phi_Bsjpsiphi_ATLAS2020B", names));
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

        meas.insert(pair<string, dato>("BRBsjpsikbst0", dato(4.13e-5, 0.12e-5, 0.07e-5, 0.14e-5, 0.45e-5))); // LHCb:2013lka

        // Angular analysis measurements: polarization fractions and CP asymmetries
        //  LHCb:2025vrr (Use only run 2 since no correlation matrix is given for the combined Run 1+2)

        TMatrixDSym CorrBsjpsikbst(7);
        CorrData.push_back(dato(0.014, 0.030)); // A_CP_0
        names.push_back("ACP_0_Bsjpsikbst0");
        CorrData.push_back(dato(-0.055, 0.066)); // A_CP_paral
        names.push_back("ACP_paral_Bsjpsikbst0");
        CorrData.push_back(dato(0.060, 0.059)); // A_CP_perp
        names.push_back("ACP_perp_Bsjpsikbst0");
        CorrData.push_back(dato(2.879, 0.087, 0. , 0., 0., true)); // delta_paral - delta_0
        names.push_back("delta_paral_Bsjpsikbst0");
        CorrData.push_back(dato(0.057, 0.068, 0. , 0., 0., true)); // delta_perp - delta_0
        names.push_back("delta_perp_Bsjpsikbst0");
        CorrData.push_back(dato(0.534, 0.015)); // f_0
        names.push_back("f_0_Bsjpsikbst0");
        CorrData.push_back(dato(0.211, 0.015)); // f_perp
        names.push_back("f_perp_Bsjpsikbst0");

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
        corrmeas.insert(pair<string, CorrelatedGaussianObservables>("Bsjpsikbst0_LHCb2025", CorrelatedGaussianObservables(CorrData, CorrBsjpsikbst)));
        corrmeas_channels.insert(pair<string, vector<string>>("Bsjpsikbst0_LHCb2025", names));
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
        //  R_Bdjpsiom_Bdjpsirho0 from LHCb:2012cw
        meas.insert(pair<string, dato>("R_Bdjpsiom_Bdjpsirho0", dato(0.86, 0.19, 0.10))); // LHCb:2012cw

        // Transversity fractions and phases from LHCb:2014vbo
        meas.insert(pair<string, dato>("f_0_Bdjpsiom", dato(0.405, 0.14, 0.035)));                                                         // LHCb:2014vbo
        meas.insert(pair<string, dato>("f_paral_Bdjpsiom", dato(0.58, 0.135, 0.035)));                                                     // LHCb:2014vbo
        meas.insert(pair<string, dato>("Delta_delta_paral_Bdjpsiom-delta_0_Bdjpsirho0", dato(123.5 / 180. * M_PI, 13.7 / 180. * M_PI, 0. , 0., 0., true)));   // LHCb:2014vbo
        meas.insert(pair<string, dato>("Delta_delta_perp_Bdjpsiom-delta_perp_Bdjpsirho0", dato(227.4 / 180. * M_PI, 84.9 / 180. * M_PI, 0., 0., 0., true))); // LHCb:2014vbo
        meas.insert(pair<string, dato>("Delta_delta_0_Bdjpsiom-delta_0_Bdjpsirho0", dato(268.8 / 180. * M_PI, 11.9 / 180. * M_PI, 0., 0., 0., true)));       // LHCb:2014vbo

        /////////////////////////////
        // Bdjpsikst0
        /////////////////////////////

        // BRBdjpsikst0
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
        pdgaverage.setName("BRBdjpsikst0");
        pdgaverage.CalculateAverage();

        meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
        data.clear();

        // Angular integrated time-dependent CP asymmetry measurements
        meas.insert(pair<string, dato>("SBdjpsikst0", dato(0.601, 0.239, 0.087))); // BaBar:2009byl
        meas.insert(pair<string, dato>("CBdjpsikst0", dato(0.025, 0.083, 0.054))); // BaBar:2009byl

        // polarization

        // f_0 Bdjpsikst0
        data.push_back(dato(0.587, 0.011));        // D0:2008nly
        data.push_back(dato(0.556, 0.009, 0.010)); // BaBar:2007rbr
        data.push_back(dato(0.562, 0.026, 0.018)); // CDF:2004dxr
        data.push_back(dato(0.574, 0.012, 0.009)); // Belle:2005qtf
        data.push_back(dato(0.59, 0.06, 0.01));    // CDF:2000edf
        data.push_back(dato(0.52, 0.07, 0.04));    // CLEO:1997ilq
        data.push_back(dato(0.65, 0.10, 0.04));    // CDF:1995kwt
        data.push_back(dato(0.97, 0.16, 0.15));    // ARGUS:1994rms

        pdgaverage.setData(data);
        pdgaverage.setName("f_0_Bdjpsikst0");
        pdgaverage.CalculateAverage();

        meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
        data.clear();

        // f_paral Bdjpsikst0
        data.push_back(dato(0.230, 0.013, 0.025)); // D0:2008nly
        data.push_back(dato(0.211, 0.010, 0.006)); // BaBar:2007rbr

        pdgaverage.setData(data);
        pdgaverage.setName("f_paral_Bdjpsikst0");
        pdgaverage.CalculateAverage();

        meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
        data.clear();

        // f_perp Bdjpsikst0
        data.push_back(dato(0.233, 0.010, 0.005)); // BaBar:2007rbr
        data.push_back(dato(0.215, 0.032, 0.006)); // CDF:2004dxr
        data.push_back(dato(0.195, 0.012, 0.008)); // Belle:2005qtf

        pdgaverage.setData(data);
        pdgaverage.setName("f_perp_Bdjpsikst0");
        pdgaverage.CalculateAverage();

        meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
        data.clear();

        // delta_paral Bdjpsikst0
        data.push_back(dato(-2.69, 0.08, 0.11, 0., 0., true)); // D0:2008nly
        data.push_back(dato(-2.93, 0.08, 0.04, 0., 0., true)); // BaBar:2007rbr

        pdgaverage.setData(data);
        pdgaverage.setName("delta_paral_Bdjpsikst0");
        pdgaverage.CalculateAverage();

        meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty(), 0., 0., 0., true)));
        data.clear();

        // delta_perp Bdjpsikst0
        data.push_back(dato(3.21, 0.06, 0.06, 0., 0., true)); // D0:2008nly
        data.push_back(dato(2.91, 0.05, 0.03, 0., 0., true)); // BaBar:2007rbr

        pdgaverage.setData(data);
        pdgaverage.setName("delta_perp_Bdjpsikst0");
        pdgaverage.CalculateAverage();

        meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty(), 0., 0., 0., true)));
        data.clear();

        // Polarization fractions and phases from LHCb:2013vga

        TMatrixDSym CorrBdjpsikst0(4);
        CorrData.push_back(dato(0.227, sqrt(0.004 * 0.004 + 0.011 * 0.011))); // f_paral
        names.push_back("f_paral_Bdjpsikst0");
        CorrData.push_back(dato(0.201, sqrt(0.004 * 0.004 + 0.008 * 0.008))); // f_perp
        names.push_back("f_perp_Bdjpsikst0");
        CorrData.push_back(dato(-2.94, sqrt(0.02 * 0.02 + 0.03 * 0.03), 0. , 0., 0., true)); // delta_paral - delta_0
        names.push_back("delta_paral_Bdjpsikst0");
        CorrData.push_back(dato(2.94, sqrt(0.02 * 0.02 + 0.02 * 0.02), 0., 0., 0., true)); // delta_perp - delta_0
        names.push_back("delta_perp_Bdjpsikst0");

        // Populate the correlation matrix
        CorrBdjpsikst0(0, 0) = 1.;
        CorrBdjpsikst0(0, 1) = CorrBdjpsikst0(1, 0) = -0.70;
        CorrBdjpsikst0(0, 2) = CorrBdjpsikst0(2, 0) = 0.12;
        CorrBdjpsikst0(0, 3) = CorrBdjpsikst0(3, 0) = 0.04;
        CorrBdjpsikst0(1, 1) = 1.;
        CorrBdjpsikst0(1, 2) = CorrBdjpsikst0(2, 1) = -0.14;
        CorrBdjpsikst0(1, 3) = CorrBdjpsikst0(3, 1) = -0.01;
        CorrBdjpsikst0(2, 2) = 1.;
        CorrBdjpsikst0(2, 3) = CorrBdjpsikst0(3, 2) = 0.64;
        CorrBdjpsikst0(3, 3) = 1.;
        // Insert correlated data into corrmeas
        corrmeas.insert(pair<string, CorrelatedGaussianObservables>("Bdjpsikst0_LHCb2013", CorrelatedGaussianObservables(CorrData, CorrBdjpsikst0)));
        corrmeas_channels.insert(pair<string, vector<string>>("Bdjpsikst0_LHCb2013", names));
        CorrData.clear();
        names.clear();

        // CP asymmetries from LHCb:2013vga
        meas.insert(pair<string, dato>("ACP_fparal_Bdjpsikst0", dato(-0.011, 0.016, 0.005)));      // LHCb:2013vga
        meas.insert(pair<string, dato>("ACP_fperp_Bdjpsikst0", dato(0.032, 0.018, 0.003)));        // LHCb:2013vga
        meas.insert(pair<string, dato>("ACP_delta_paral_Bdjpsikst0", dato(-0.003, 0.007, 0.002))); // LHCb:2013vga
        meas.insert(pair<string, dato>("ACP_delta_perp_Bdjpsikst0", dato(0.003, 0.005, 0.001)));   // LHCb:2013vga

        // CBdjpsikst0
        meas.insert(pair<string, dato>("CBdjpsikst0", dato(0.025, 0.083, 0.054))); // BaBar:2009byl

        // SBdjpsikst0
        meas.insert(pair<string, dato>("SBdjpsikst0", dato(0.601, 0.239, 0.087))); // BaBar:2009byl

        /////////////////////////////
        // Bdjpsirho0
        /////////////////////////////

        // BRBdjpsirho0
        data.push_back(dato(2.515e-5, 0.10e-5, 0.165e-5)); // LHCb:2014vbo
        data.push_back(dato(2.7e-5, 0.3e-5, 0.2e-5));      // BaBar:2007yvx

        pdgaverage.setData(data);
        pdgaverage.setName("BRBdjpsirho0");
        pdgaverage.CalculateAverage();

        meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
        data.clear();

        // CP violation for the different polarizations of Bdjpsirho0 from LHCb:2014xpr

        CorrData.push_back(dato(-0.047, 0.034, 0.011)); // Alpha_0
        names.push_back("alpha_CP_0_Bdjpsirho0");
        CorrData.push_back(dato(-0.060, 0.060, 0.007)); // Alpha_paral
        names.push_back("alpha_CP_paral_Bdjpsirho0");
        CorrData.push_back(dato(0.020, 0.109, 0.018)); // Alpha_perp
        names.push_back("alpha_CP_perp_Bdjpsirho0");
        CorrData.push_back(dato(-3.3 / 180. * M_PI, 7.2 / 180. * M_PI, 1.7 / 180. * M_PI)); // 2beta_perp - 2beta_0
        names.push_back("2beta_perp_Bdjpsirho0");
        CorrData.push_back(dato(-0.5 / 180. * M_PI, 6.5 / 180. * M_PI, 1.6 / 180. * M_PI)); // 2beta_paral - 2beta_0
        names.push_back("2beta_paral_Bdjpsirho0");
        CorrData.push_back(dato(42.1 / 180. * M_PI, 10.2 / 180. * M_PI, 5.0 / 180. * M_PI)); // 2beta_0
        names.push_back("2beta_0_Bdjpsirho0");
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
        corrmeas_channels.insert(pair<string, vector<string>>("Bdjpsirh_LHCb2014xpr", names));
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

        // BRBpjpsikstp
        meas.insert(pair<string, dato>("BRBpjpsikstp", dato(1.43e-3, 0.08e-3))); // PDGlive 1/2025

        // polarization fractions from Belle:2005qtf
        meas.insert(pair<string, dato>("f_0_Bpjpsikstp", dato(0.604, 0.015, 0.018)));    // Belle:2005qtf
        meas.insert(pair<string, dato>("f_perp_Bpjpsikstp", dato(0.180, 0.014, 0.010))); // Belle:2005qtf
        // CP asymmetry from BaBar:2004htr
        meas.insert(pair<string, dato>("ACPBpjpsikstp", dato(0.048, 0.029, 0.016))); // BaBar:2004htr

        /////////////////////////////
        // Bpjpsirhop
        /////////////////////////////

        // BRBpjpsirhop
        data.push_back(dato(3.81e-5, 0.25e-5, 0.35e-5)); // LHCb:2018pil
        data.push_back(dato(5.0e-5, 0.7e-5, 0.31e-5));   // BaBar:2007yvx

        pdgaverage.setData(data);
        pdgaverage.setName("BRBpjpsirhop");
        pdgaverage.CalculateAverage();
        meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
        data.clear();

        // CP asymmetry

        data.push_back(dato(-0.045, 0.056, 0.008)); // LHCb:2018pil
        data.push_back(dato(-0.11, 0.12, 0.08));    // BaBar:2007yvx

        pdgaverage.setData(data);
        pdgaverage.setName("ACPBpjpsirhop");
        pdgaverage.CalculateAverage();
        meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
        data.clear();
    }
    if (BDDb)
    {

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

        // BRBsdspdsm (Bs→Ds⁺Ds⁻)
        data.push_back(dato(4.0e-3, 0.2e-3, 0.5e-3)); // LHCb:2013sad
        data.push_back(dato(5.9e-3, 1.0e-3, 1.3e-3)); // Belle:2012tsw
        data.push_back(dato(5.4e-3, 0.8e-3, 0.8e-3)); // CDF:2012xmd

        pdgaverage.setData(data);
        pdgaverage.setName("BRBsdspdsm");
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

        // BRBpdspd0b (B⁺→D⁺D̄⁰ₛ)
        data.push_back(dato(8.6e-3, 0.2e-3, 1.1e-3)); // LHCb:2013sad
        data.push_back(dato(9.5e-3, 2.0e-3, 0.8e-3)); // BaBar:2006jvx
        data.push_back(dato(9.8e-3, 2.6e-3, 0.9e-3)); // CLEO:1995psi
        data.push_back(dato(14e-3, 8e-3, 1e-3));      // ARGUS:1991xej
        data.push_back(dato(13e-3, 6e-3, 1e-3));      // CLEO:1990mqz

        pdgaverage.setData(data);
        pdgaverage.setName("BRBpdspd0b");
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
        corrmeas_channels.insert(pair<string, vector<string>>("CS_Bddpdm_LHCb2024", names));
        names.clear();
        CorrData.clear();

        // Bsdspdsm : C and S observables from LHCb Run2 (part of arXiv:2409.03009 combination)
        // We use C and S instead of phi_s and lambda to avoid using pre-averaged values
        // From LHCb Run2: C = 0.128 ± 0.103(stat) ± 0.010(syst), S = 0.552 ± 0.100(stat) ± 0.010(syst)
        // NOTE: These are for Bs→Ds⁺Ds⁻ from the same paper as Bd→D⁺D⁻
        CorrData.push_back(dato(-0.053, 0.096, 0.020)); // C observable (NOTE: Check paper for actual values)
        names.push_back("C_Bsdspdsm");
        CorrData.push_back(dato(-0.055, 0.092, 0.021)); // S observable (NOTE: Check paper for actual values)
        names.push_back("S_Bsdspdsm");

        CorrStat(0, 0) = 1.;
        CorrStat(1, 1) = 1.;
        CorrStat(0, 1) = 0.; // Correlation coefficient (CHECK PAPER)
        CorrStat(1, 0) = 0.;

        CorrSyst(0, 0) = 1.;
        CorrSyst(1, 1) = 1.;
        CorrSyst(0, 1) = 0.;
        CorrSyst(1, 0) = 0.;

        corrmeas.insert(pair<string, CorrelatedGaussianObservables>("CS_Bsdspdsm_LHCb2024", CorrelatedGaussianObservables(CorrData, CorrStat, CorrSyst)));
        corrmeas_channels.insert(pair<string, vector<string>>("CS_Bsdspdsm_LHCb2024", names));
        names.clear();
        CorrData.clear();

        // ACPBpdpd0b and ACPBpdspd0b from LHCb:2023wbb (correlated)
        dato Bpdpd0b_acp(2.5e-2, 1.0e-2, 0.4e-2, 0.3e-2);
        names.push_back("ACPBpdpd0b");
        dato Bpdspd0b_acp(0.5e-2, 0.2e-2, 0.5e-2, 0.3e-2);
        names.push_back("ACPBpdspd0b");

        CorrData.push_back(dato(Bpdpd0b_acp.getMean(), Bpdpd0b_acp.getSigma()));   // ACPBpdpd0b
        CorrData.push_back(dato(Bpdspd0b_acp.getMean(), Bpdspd0b_acp.getSigma())); // ACPBpdspd0b

        TMatrixDSym Corr_ACP_charged(2);
        Corr_ACP_charged(0, 0) = 1.0;
        Corr_ACP_charged(1, 1) = 1.0;
        Corr_ACP_charged(0, 1) = 0.386; // Correlation from LHCb:2023wbb
        Corr_ACP_charged(1, 0) = 0.386;

        corrmeas.insert(pair<string, CorrelatedGaussianObservables>("ACP_Bpdpd0b_Bpdspd0b_LHCb2023", CorrelatedGaussianObservables(CorrData, Corr_ACP_charged)));
        corrmeas_channels.insert(pair<string, vector<string>>("ACP_Bpdpd0b_Bpdspd0b_LHCb2023", names));
        names.clear();
        CorrData.clear();
    }
    cout << "All meas inserted" << endl;

    // Create histograms for all observables

    vector<string> paramsfor2dhistos;

    for (const auto &channel : channelNamesSU3)
    {
        // create histos for mod and phase as well as real and imaginary parts for the effective parameters
        for (const auto &param : channelParameters[channel])
        {
            histos.createH1D(param, 500, 0.0, 0.0);
            string newStr;
            size_t length = param.length();

            newStr.reserve(length + 1);

            if (length >= 3 && param.substr(length - 3) == "_re")
            {
                newStr.append(param, 0, length - 3); // Append text before "_re"
                // strip the possible trailing "delta_" from the parameter name to get the full parameter
                newStr.erase(0, newStr.find("delta_") == 0 ? 6 : 0); // Remove "delta_" if it exists at the start
                histos.createH1D(newStr + "_abs", 500, 0.0, 0.0);
                paramsfor2dhistos.push_back(newStr + "_abs");
            }
            else if (length >= 3 && param.substr(length - 3) == "_im")
            {
                newStr.append(param, 0, length - 3);                 // Append text before "_im"
                newStr.erase(0, newStr.find("delta_") == 0 ? 6 : 0); // Remove "delta_" if it exists at the start
                histos.createH1D(newStr + "_arg", 500, 0.0, 0.0);
                paramsfor2dhistos.push_back(newStr + "_arg");
            }
        }
    }

    for (size_t i = 0; i < paramsfor2dhistos.size(); ++i)
        for (size_t j = i + 1; j < paramsfor2dhistos.size(); ++j)
        {
            const auto &param = paramsfor2dhistos[i];
            const auto &param2 = paramsfor2dhistos[j];
            histos.createH2D(param, param2, 500, 0.0, 0.0, 500, 0.0, 0.0);
        }

    for (const auto &channel : channels)
    {

        histos.createH1D("BR_" + channel, 500, 0.0, 0.0);

        auto parsedChannel = parseChannel(channel);

        if (parsedChannel.first == "Bp")
        {
            histos.createH1D("ACP_" + channel, 500, 0.0, 0.0);
        }
        else
        {

            histos.createH1D("C_" + channel, 500, 0.0, 0.0);

            histos.createH1D("S_" + channel, 500, 0.0, 0.0);
            histos.createH1D("DeltaS_" + channel, 500, 0.0, 0.0);
        }
    }

    for (const auto &corr_entry : corrmeas_channels)
    {
        const string &corr_name = corr_entry.first;
        const vector<string> &obs_names = corr_entry.second;

        for (const auto &obs_name : obs_names)
        {
            histos.createH1D(obs_name, 500, 0.0, 0.0);
            if (obs_name.find("S_") != string::npos || obs_name.find("phis") != string::npos || obs_name.find("2beta") != string::npos)
            {
                histos.createH1D("Delta" + obs_name, 500, 0.0, 0.0);
            }
        }
    }

    histos.createH1D("2beta", 500, 0.0, 0.0);
    histos.createH1D("phis", 500, 0.0, 0.0);
    histos.createH1D("myphis", 500, 0.0, 0.0);
    histos.createH1D("myphid", 500, 0.0, 0.0);

    histos.createH1D("LogLikelihood", 500, 0.0, 0.0);
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
        vector<string> params = {
            "E2t_ccsd_BJPSIP_re",
            "E2t_ccsd_BJPSIP_im",
            "R_G2t_scd_BJPSIP_re",
            "R_G2t_scd_BJPSIP_im",
        };
        channelParameters[channel] = params;

        addAmplitudeParameter("E2t_ccsd_BJPSIP_re", 0., 5.);
        addAmplitudeParameter("E2t_ccsd_BJPSIP_im", 0., 0.);
        addAmplitudeParameter("R_G2t_scd_BJPSIP_re", -1., 1.1);
        addAmplitudeParameter("R_G2t_scd_BJPSIP_im", -1., 1.1);
    }
    else if (channel == "Bdjpsip0")
    {
        vector<string> params = {
            "delta_E2t_ccdd_BJPSIP_re",
            "delta_E2t_ccdd_BJPSIP_im",
            "dP4EW_ucd_BPJPSI_re",
            "dP4EW_ucd_BPJPSI_im",
            "R_EA2_ddcd_BPJPSI_re",
            "R_EA2_ddcd_BPJPSI_im",
            "delta_G2t_dcd_BJPSIP_re",
            "delta_G2t_dcd_BJPSIP_im"};
        channelParameters[channel] = params;

        addAmplitudeParameter("dP4EW_ucd_BPJPSI_re", -ewp_limit, ewp_limit);
        addAmplitudeParameter("dP4EW_ucd_BPJPSI_im", -ewp_limit, ewp_limit);
        addAmplitudeParameter("R_EA2_ddcd_BPJPSI_re", -1., 1.1);
        addAmplitudeParameter("R_EA2_ddcd_BPJPSI_im", -1., 1.1);
        addSU3BreakingParameter("delta_E2t_ccdd_BJPSIP_re", "E2t_ccsd_BJPSIP_re");
        addSU3BreakingParameter("delta_E2t_ccdd_BJPSIP_im", "E2t_ccsd_BJPSIP_im");
        addSU3BreakingParameter("delta_G2t_dcd_BJPSIP_re", "G2t_scd_BJPSIP_re");
        addSU3BreakingParameter("delta_G2t_dcd_BJPSIP_im", "G2t_scd_BJPSIP_im");
    }
    else if (channel == "Bdjpsieta8")
    {
        vector<string> params = {
            "delta_E2t_ccdd_BJPSIP_re",
            "delta_E2t_ccdd_BJPSIP_im",
            "dP4EW_ucd_BPJPSI_re",
            "dP4EW_ucd_BPJPSI_im",
            "R_EA2t_ccdd_BJPSIP_re",
            "R_EA2t_ccdd_BJPSIP_im",
            "delta_EA2t_ccsd_BJPSIP_re",
            "delta_EA2t_ccsd_BJPSIP_im",
            "R_EA2_ddcd_BPJPSI_re",
            "R_EA2_ddcd_BPJPSI_im",
            "delta_G2t_dcd_BJPSIP_re",
            "delta_G2t_dcd_BJPSIP_im",
            "R_G4t_cdd_BJPSIP_re",
            "R_G4t_cdd_BJPSIP_im",
            "delta_G4t_csd_BJPSIP_re",
            "delta_G4t_csd_BJPSIP_im"};
        channelParameters[channel] = params;

        addSU3BreakingParameter("delta_E2t_ccdd_BJPSIP_re", "E2t_ccsd_BJPSIP_re");
        addSU3BreakingParameter("delta_E2t_ccdd_BJPSIP_im", "E2t_ccsd_BJPSIP_im");
        addAmplitudeParameter("dP4EW_ucd_BPJPSI_re", -ewp_limit, ewp_limit);
        addAmplitudeParameter("dP4EW_ucd_BPJPSI_im", -ewp_limit, ewp_limit);
        addAmplitudeParameter("R_EA2t_ccdd_BJPSIP_re", -1., 1.1);
        addAmplitudeParameter("R_EA2t_ccdd_BJPSIP_im", -1., 1.1);
        addSU3BreakingParameter("delta_EA2t_ccsd_BJPSIP_re", "EA2t_ccdd_BJPSIP_re");
        addSU3BreakingParameter("delta_EA2t_ccsd_BJPSIP_im", "EA2t_ccdd_BJPSIP_im");
        addAmplitudeParameter("R_EA2_ddcd_BPJPSI_re", -1., 1.1);
        addAmplitudeParameter("R_EA2_ddcd_BPJPSI_im", -1., 1.1);
        addSU3BreakingParameter("delta_G2t_dcd_BJPSIP_re", "G2t_scd_BJPSIP_re");
        addSU3BreakingParameter("delta_G2t_dcd_BJPSIP_im", "G2t_scd_BJPSIP_im");
        addAmplitudeParameter("R_G4t_cdd_BJPSIP_re", -1., 1.1);
        addAmplitudeParameter("R_G4t_cdd_BJPSIP_im", -1., 1.1);
        addSU3BreakingParameter("delta_G4t_csd_BJPSIP_re", "G4t_cdd_BJPSIP_re");
        addSU3BreakingParameter("delta_G4t_csd_BJPSIP_im", "G4t_cdd_BJPSIP_im");
    }
    else if (channel == "Bdjpsieta1")
    {
        vector<string> params = {
            "delta_E2t_ccdd_BJPSIP_re",
            "delta_E2t_ccdd_BJPSIP_im",
            "dP4EW_ucd_BPJPSI_re",
            "dP4EW_ucd_BPJPSI_im",
            "R_EA2t_ccdd_BJPSIP_re",
            "R_EA2t_ccdd_BJPSIP_im",
            "delta_EA2t_ccsd_BJPSIP_re",
            "delta_EA2t_ccsd_BJPSIP_im",
            "R_EA2_ddcd_BPJPSI_re",
            "R_EA2_ddcd_BPJPSI_im",
            "delta_G2t_dcd_BJPSIP_re",
            "delta_G2t_dcd_BJPSIP_im",
            "R_G4t_cdd_BJPSIP_re",
            "R_G4t_cdd_BJPSIP_im",
            "delta_G4t_csd_BJPSIP_re",
            "delta_G4t_csd_BJPSIP_im"};
        channelParameters[channel] = params;

        addSU3BreakingParameter("delta_E2t_ccdd_BJPSIP_re", "E2t_ccsd_BJPSIP_re");
        addSU3BreakingParameter("delta_E2t_ccdd_BJPSIP_im", "E2t_ccsd_BJPSIP_im");
        addAmplitudeParameter("dP4EW_ucd_BPJPSI_re", -ewp_limit, ewp_limit);
        addAmplitudeParameter("dP4EW_ucd_BPJPSI_im", -ewp_limit, ewp_limit);
        addAmplitudeParameter("R_EA2t_ccdd_BJPSIP_re", -1., 1.1);
        addAmplitudeParameter("R_EA2t_ccdd_BJPSIP_im", -1., 1.1);
        addSU3BreakingParameter("delta_EA2t_ccsd_BJPSIP_re", "EA2t_ccdd_BJPSIP_re");
        addSU3BreakingParameter("delta_EA2t_ccsd_BJPSIP_im", "EA2t_ccdd_BJPSIP_im");
        addAmplitudeParameter("R_EA2_ddcd_BPJPSI_re", -1., 1.1);
        addAmplitudeParameter("R_EA2_ddcd_BPJPSI_im", -1., 1.1);
        addSU3BreakingParameter("delta_G2t_dcd_BJPSIP_re", "G2t_scd_BJPSIP_re");
        addSU3BreakingParameter("delta_G2t_dcd_BJPSIP_im", "G2t_scd_BJPSIP_im");
        addAmplitudeParameter("R_G4t_cdd_BJPSIP_re", -1., 1.1);
        addAmplitudeParameter("R_G4t_cdd_BJPSIP_im", -1., 1.1);
        addSU3BreakingParameter("delta_G4t_csd_BJPSIP_re", "G4t_cdd_BJPSIP_re");
        addSU3BreakingParameter("delta_G4t_csd_BJPSIP_im", "G4t_cdd_BJPSIP_im");
    }
    else if (channel == "Bpjpsikp")
    {
        vector<string> params = {
            "E2t_ccsd_BJPSIP_re",
            "E2t_ccsd_BJPSIP_im",
            "R_G2t_scd_BJPSIP_re",
            "R_G2t_scd_BJPSIP_im",
            "dP2EW_scu_BPJPSI_re",
            "dP2EW_scu_BPJPSI_im",
            "R_EA1_sdcd_BPJPSI_re",
            "R_EA1_sdcd_BPJPSI_im"};
        channelParameters[channel] = params;
        addAmplitudeParameter("E2t_ccsd_BJPSIP_re", 0., 10.);
        addAmplitudeParameter("E2t_ccsd_BJPSIP_im", 0., 0.);
        addAmplitudeParameter("R_G2t_scd_BJPSIP_re", -1., 1.1);
        addAmplitudeParameter("R_G2t_scd_BJPSIP_im", -1., 1.1);
        addAmplitudeParameter("dP2EW_scu_BPJPSI_re", -ewp_limit, ewp_limit);
        addAmplitudeParameter("dP2EW_scu_BPJPSI_im", -ewp_limit, ewp_limit);
        addAmplitudeParameter("R_EA1_sdcd_BPJPSI_re", -1., 1.1);
        addAmplitudeParameter("R_EA1_sdcd_BPJPSI_im", -1., 1.1);
    }
    else if (channel == "Bpjpsipp")
    {
        vector<string> params = {
            "delta_E2t_ccdd_BJPSIP_re",
            "delta_E2t_ccdd_BJPSIP_im",
            "delta_G2t_dcd_BJPSIP_re",
            "delta_G2t_dcd_BJPSIP_im",
            "delta_dP2EW_dcu_BPJPSI_re",
            "delta_dP2EW_dcu_BPJPSI_im",
            "delta_EA1_ddcd_BPJPSI_re",
            "delta_EA1_ddcd_BPJPSI_im"};
        channelParameters[channel] = params;
        addSU3BreakingParameter("delta_E2t_ccdd_BJPSIP_re", "E2t_ccsd_BJPSIP_re");
        addSU3BreakingParameter("delta_E2t_ccdd_BJPSIP_im", "E2t_ccsd_BJPSIP_im");
        addSU3BreakingParameter("delta_G2t_dcd_BJPSIP_re", "G2t_scd_BJPSIP_re");
        addSU3BreakingParameter("delta_G2t_dcd_BJPSIP_im", "G2t_scd_BJPSIP_im");
        addSU3BreakingParameter("delta_dP2EW_dcu_BPJPSI_re", "dP2EW_scu_BPJPSI_re");
        addSU3BreakingParameter("delta_dP2EW_dcu_BPJPSI_im", "dP2EW_scu_BPJPSI_im");
        addSU3BreakingParameter("delta_EA1_ddcd_BPJPSI_re", "EA1_sdcd_BPJPSI_re");
        addSU3BreakingParameter("delta_EA1_ddcd_BPJPSI_im", "EA1_sdcd_BPJPSI_im");
    }
    else if (channel == "Bsjpsip0")
    {
        vector<string> params = {
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
        vector<string> params = {
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
        addSU3BreakingParameter("delta_E2t_ccds_BJPSIP_re", "E2t_ccdd_BJPSIP_re");
        addSU3BreakingParameter("delta_E2t_ccds_BJPSIP_im", "E2t_ccdd_BJPSIP_im");
        addSU3BreakingParameter("delta_G2t_dcd_BJPSIP_re", "G2t_scd_BJPSIP_re");
        addSU3BreakingParameter("delta_G2t_dcd_BJPSIP_im", "G2t_scd_BJPSIP_im");
        addSU3BreakingParameter("delta_G2t_dcs_BJPSIP_re", "G2t_dcd_BJPSIP_re");
        addSU3BreakingParameter("delta_G2t_dcs_BJPSIP_im", "G2t_dcd_BJPSIP_im");
    }
    else if (channel == "Bsjpsieta8")
    {
        vector<string> params = {
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
        addSU3BreakingParameter("delta_G4t_css_BJPSIP_re", "G4t_cds_BJPSIP_re");
        addSU3BreakingParameter("delta_G4t_css_BJPSIP_im", "G4t_cds_BJPSIP_im");
        addSU3BreakingParameter("delta_EA2t_ccds_BJPSIP_re", "EA2t_ccdd_BJPSIP_re");
        addSU3BreakingParameter("delta_EA2t_ccds_BJPSIP_im", "EA2t_ccdd_BJPSIP_im");
        addSU3BreakingParameter("delta_EA2t_ccss_BJPSIP_re", "EA2t_ccds_BJPSIP_re");
        addSU3BreakingParameter("delta_EA2t_ccss_BJPSIP_im", "EA2t_ccds_BJPSIP_im");
    }
    else if (channel == "Bsjpsieta1")
    {
        vector<string> params = {
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
        addSU3BreakingParameter("delta_G4t_css_BJPSIP_re", "G4t_cds_BJPSIP_re");
        addSU3BreakingParameter("delta_G4t_css_BJPSIP_im", "G4t_cds_BJPSIP_im");
        addSU3BreakingParameter("delta_EA2t_ccds_BJPSIP_re", "EA2t_ccdd_BJPSIP_re");
        addSU3BreakingParameter("delta_EA2t_ccds_BJPSIP_im", "EA2t_ccdd_BJPSIP_im");
        addSU3BreakingParameter("delta_EA2t_ccss_BJPSIP_re", "EA2t_ccds_BJPSIP_re");
        addSU3BreakingParameter("delta_EA2t_ccss_BJPSIP_im", "EA2t_ccds_BJPSIP_im");
    }
    // Vector channels with helicity amplitudes
    else if (channel == "Bsjpsiphi")
    {
        vector<string> params = {
            "E2t_ccss_BJPSIV_0_re",
            "E2t_ccss_BJPSIV_0_im",
            "R_G2t_scs_BJPSIV_0_re",
            "R_G2t_scs_BJPSIV_0_im",
            "R_EA2t_ccss_BJPSIV_0_re",
            "R_EA2t_ccss_BJPSIV_0_im",
            "R_G4t_css_BJPSIV_0_re",
            "R_G4t_css_BJPSIV_0_im",
            "E2t_ccss_BJPSIV_paral_re",
            "E2t_ccss_BJPSIV_paral_im",
            "R_G2t_scs_BJPSIV_paral_re",
            "R_G2t_scs_BJPSIV_paral_im",
            "R_EA2t_ccss_BJPSIV_paral_re",
            "R_EA2t_ccss_BJPSIV_paral_im",
            "R_G4t_css_BJPSIV_paral_re",
            "R_G4t_css_BJPSIV_paral_im",
            "E2t_ccss_BJPSIV_perp_re",
            "E2t_ccss_BJPSIV_perp_im",
            "R_G2t_scs_BJPSIV_perp_re",
            "R_G2t_scs_BJPSIV_perp_im",
            "R_EA2t_ccss_BJPSIV_perp_re",
            "R_EA2t_ccss_BJPSIV_perp_im",
            "R_G4t_css_BJPSIV_perp_re",
            "R_G4t_css_BJPSIV_perp_im"};
        channelParameters[channel] = params;

        addAmplitudeParameter("E2t_ccss_BJPSIV_re", 0., 10., true);
        addAmplitudeParameter("E2t_ccss_BJPSIV_im", 0., 0., true);
        addAmplitudeParameter("R_G2t_scs_BJPSIV_re", -1., 1.1, true);
        addAmplitudeParameter("R_G2t_scs_BJPSIV_im", -1., 1.1, true);
        addAmplitudeParameter("R_EA2t_ccss_BJPSIV_re", -1., 1.1, true);
        addAmplitudeParameter("R_EA2t_ccss_BJPSIV_im", -1., 1.1, true);
        addAmplitudeParameter("R_G4t_css_BJPSIV_re", -1., 1.1, true);
        addAmplitudeParameter("R_G4t_css_BJPSIV_im", -1., 1.1, true);
    }
    else if (channel == "Bsjpsiom")
    {
        vector<string> params = {
            "delta_EA2t_ccds_BJPSIV_0_re",
            "delta_EA2t_ccds_BJPSIV_0_im",
            "delta_G4t_cds_BJPSIV_0_re",
            "delta_G4t_cds_BJPSIV_0_im",
            "dP4EW_ucs_BVJPSI_0_re",
            "dP4EW_ucs_BVJPSI_0_im",
            "R_EA2_ddcs_BVJPSI_0_re",
            "R_EA2_ddcs_BVJPSI_0_im",
            "delta_EA2t_ccds_BJPSIV_paral_re",
            "delta_EA2t_ccds_BJPSIV_paral_im",
            "delta_G4t_cds_BJPSIV_paral_re",
            "delta_G4t_cds_BJPSIV_paral_im",
            "dP4EW_ucs_BVJPSI_paral_re",
            "dP4EW_ucs_BVJPSI_paral_im",
            "R_EA2_ddcs_BVJPSI_paral_re",
            "R_EA2_ddcs_BVJPSI_paral_im",
            "delta_EA2t_ccds_BJPSIV_perp_re",
            "delta_EA2t_ccds_BJPSIV_perp_im",
            "delta_G4t_cds_BJPSIV_perp_re",
            "delta_G4t_cds_BJPSIV_perp_im",
            "dP4EW_ucs_BVJPSI_perp_re",
            "dP4EW_ucs_BVJPSI_perp_im",
            "R_EA2_ddcs_BVJPSI_perp_re",
            "R_EA2_ddcs_BVJPSI_perp_im"};
        channelParameters[channel] = params;
        addSU3BreakingParameter("delta_EA2t_ccds_BJPSIV_re", "EA2t_ccss_BJPSIV_re", true);
        addSU3BreakingParameter("delta_EA2t_ccds_BJPSIV_im", "EA2t_ccss_BJPSIV_im", true);
        addSU3BreakingParameter("delta_G4t_cds_BJPSIV_re", "G4t_css_BJPSIV_re", true);
        addSU3BreakingParameter("delta_G4t_cds_BJPSIV_im", "G4t_css_BJPSIV_im", true);
        addAmplitudeParameter("dP4EW_ucs_BVJPSI_re", -ewp_limit, ewp_limit, true);
        addAmplitudeParameter("dP4EW_ucs_BVJPSI_im", -ewp_limit, ewp_limit, true);
        addAmplitudeParameter("R_EA2_ddcs_BVJPSI_re", -1., 1.1, true);
        addAmplitudeParameter("R_EA2_ddcs_BVJPSI_im", -1., 1.1, true);
    }
    else if (channel == "Bsjpsikbst0")
    {
        vector<string> params = {
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
        vector<string> params = {
            "dP4EW_ucs_BVJPSI_0_re",
            "dP4EW_ucs_BVJPSI_0_im",
            "R_EA2_ddcs_BVJPSI_0_re",
            "R_EA2_ddcs_BVJPSI_0_im",
            "dP4EW_ucs_BVJPSI_paral_re",
            "dP4EW_ucs_BVJPSI_paral_im",
            "R_EA2_ddcs_BVJPSI_paral_re",
            "R_EA2_ddcs_BVJPSI_paral_im",
            "dP4EW_ucs_BVJPSI_perp_re",
            "dP4EW_ucs_BVJPSI_perp_im",
            "R_EA2_ddcs_BVJPSI_perp_re",
            "R_EA2_ddcs_BVJPSI_perp_im"};
        channelParameters[channel] = params;

        addAmplitudeParameter("dP4EW_ucs_BVJPSI_re", -ewp_limit, ewp_limit, true);
        addAmplitudeParameter("dP4EW_ucs_BVJPSI_im", -ewp_limit, ewp_limit, true);
        addAmplitudeParameter("R_EA2_ddcs_BVJPSI_re", -1., 1.1, true);
        addAmplitudeParameter("R_EA2_ddcs_BVJPSI_im", -1., 1.1, true);
    }
    else if (channel == "Bdjpsiom")
    {
        vector<string> params = {
            "delta_E2t_ccds_BJPSIV_0_re",
            "delta_E2t_ccds_BJPSIV_0_im",
            "delta_E2t_ccdd_BJPSIV_0_re",
            "delta_E2t_ccdd_BJPSIV_0_im",
            "delta_G2t_dcs_BJPSIV_0_re",
            "delta_G2t_dcs_BJPSIV_0_im",
            "delta_G2t_dcd_BJPSIV_0_re",
            "delta_G2t_dcd_BJPSIV_0_im",
            "delta_dP4EW_ucd_BVJPSI_0_re",
            "delta_dP4EW_ucd_BVJPSI_0_im",
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
            "delta_dP4EW_ucd_BVJPSI_paral_re",
            "delta_dP4EW_ucd_BVJPSI_paral_im",
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
            "delta_dP4EW_ucd_BVJPSI_perp_re",
            "delta_dP4EW_ucd_BVJPSI_perp_im",
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
        addSU3BreakingParameter("delta_dP4EW_ucd_BVJPSI_re", "dP4EW_ucs_BVJPSI_re", true);
        addSU3BreakingParameter("delta_dP4EW_ucd_BVJPSI_im", "dP4EW_ucs_BVJPSI_im", true);
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
    else if (channel == "Bdjpsikst0")
    {
        vector<string> params = {
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
        vector<string> params = {
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
    else if (channel == "Bdjpsiphi")
    {
        vector<string> params = {
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
    else if (channel == "Bpjpsikstp")
    {
        vector<string> params = {
            "delta_E2t_ccsd_BJPSIV_0_re",
            "delta_E2t_ccsd_BJPSIV_0_im",
            "dP2EW_scu_BJPSIV_0_re",
            "dP2EW_scu_BJPSIV_0_im",
            "R_EA1_sdcd_BVJPSI_0_re",
            "R_EA1_sdcd_BVJPSI_0_im",
            "delta_G2t_scd_BJPSIV_0_re",
            "delta_G2t_scd_BJPSIV_0_im",
            "delta_E2t_ccsd_BJPSIV_paral_re",
            "delta_E2t_ccsd_BJPSIV_paral_im",
            "dP2EW_scu_BJPSIV_paral_re",
            "dP2EW_scu_BJPSIV_paral_im",
            "R_EA1_sdcd_BVJPSI_paral_re",
            "R_EA1_sdcd_BVJPSI_paral_im",
            "delta_G2t_scd_BJPSIV_paral_re",
            "delta_G2t_scd_BJPSIV_paral_im",
            "delta_E2t_ccsd_BJPSIV_perp_re",
            "delta_E2t_ccsd_BJPSIV_perp_im",
            "dP2EW_scu_BJPSIV_perp_re",
            "dP2EW_scu_BJPSIV_perp_im",
            "R_EA1_sdcd_BVJPSI_perp_re",
            "R_EA1_sdcd_BVJPSI_perp_im",
            "delta_G2t_scd_BJPSIV_perp_re",
            "delta_G2t_scd_BJPSIV_perp_im"};

        channelParameters[channel] = params;

        addSU3BreakingParameter("delta_E2t_ccsd_BJPSIV_re", "E2t_ccss_BJPSIV_re", true);
        addSU3BreakingParameter("delta_E2t_ccsd_BJPSIV_im", "E2t_ccss_BJPSIV_im", true);
        addAmplitudeParameter("dP2EW_scu_BJPSIV_re", -ewp_limit, ewp_limit, true);
        addAmplitudeParameter("dP2EW_scu_BJPSIV_im", -ewp_limit, ewp_limit, true);
        addAmplitudeParameter("R_EA1_sdcd_BVJPSI_re", -10., 11., true);
        addAmplitudeParameter("R_EA1_sdcd_BVJPSI_im", -10., 11., true);
        addSU3BreakingParameter("delta_G2t_scd_BJPSIV_re", "G2t_scs_BJPSIV_re", true);
        addSU3BreakingParameter("delta_G2t_scd_BJPSIV_im", "G2t_scs_BJPSIV_im", true);
    }
    else if (channel == "Bpjpsirhop")
    {
        vector<string> params = {
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
    else if (channel == "Bsdspdsm")
    {
        // Bs -> Ds+ Ds- parameters (base channel for s->c)
        // b → c(c̄s), spectator s
        vector<string> params = {
            "E1t_sccs_BDDb_re", "E1t_sccs_BDDb_im",
            "R_A2t_cscs_BDbD_re", "R_A2t_cscs_BDbD_im",
            "R_G1t_scs_BDDb_re", "R_G1t_scs_BDDb_im",
            "R_G3t_css_BDDb_re", "R_G3t_css_BDDb_im"};
        channelParameters[channel] = params;

        addAmplitudeParameter("E1t_sccs_BDDb_re", 0., 20.);
        addAmplitudeParameter("E1t_sccs_BDDb_im", 0., 0.);
        addAmplitudeParameter("R_A2t_cscs_BDbD_re", -1., 1.1);
        addAmplitudeParameter("R_A2t_cscs_BDbD_im", -1., 1.1);
        addAmplitudeParameter("R_G1t_scs_BDDb_re", -1., 1.1);
        addAmplitudeParameter("R_G1t_scs_BDDb_im", -1., 1.1);
        addAmplitudeParameter("R_G3t_css_BDDb_re", -1., 1.1);
        addAmplitudeParameter("R_G3t_css_BDDb_im", -1., 1.1);
    }
    else if (channel == "Bsdpdsm")
    {
        // b → c(c̄d), spectator s
        vector<string> params = {
            "delta_E1t_dccs_BDDb_re", "delta_E1t_dccs_BDDb_im",
            "delta_G1t_dcs_BDDb_re", "delta_G1t_dcs_BDDb_im"};
        channelParameters[channel] = params;

        addSU3BreakingParameter("delta_E1t_dccs_BDDb_re", "E1t_sccs_BDDb_re");
        addSU3BreakingParameter("delta_E1t_dccs_BDDb_im", "E1t_sccs_BDDb_im");
        addSU3BreakingParameter("delta_G1t_dcs_BDDb_re", "G1t_scs_BDDb_re");
        addSU3BreakingParameter("delta_G1t_dcs_BDDb_im", "G1t_scs_BDDb_im");
    }
    else if (channel == "Bsdpdm")
    {
        // b → c(c̄s), spectator s
        vector<string> params = {
            "delta_A2t_cdcs_BDbD_re", "delta_A2t_cdcs_BDbD_im",
            "delta_G3_cds_BDDb_re", "delta_G3_cds_BDDb_im"};
        channelParameters[channel] = params;

        addSU3BreakingParameter("delta_A2t_cdcs_BDbD_re", "A2t_cscs_BDbD_re");
        addSU3BreakingParameter("delta_A2t_cdcs_BDbD_im", "A2t_cscs_BDbD_im");
        addSU3BreakingParameter("delta_G3_cds_BDDb_re", "G3t_css_BDDb_re");
        addSU3BreakingParameter("delta_G3_cds_BDDb_im", "G3t_css_BDDb_im");
    }
    else if (channel == "Bsd0d0b")
    {
        // b → c(c̄s), spectator s
        vector<string> params = {
            "delta_A2t_cdcs_BDbD_re", "delta_A2t_cdcs_BDbD_im",
            "dP3EW_ucs_BDbD_re", "dP3EW_ucs_BDbD_im",
            "R_A2_dcds_BDDb_re", "R_A2_dcds_BDDb_im",
            "delta_G3t_cds_BDDb_re", "delta_G3t_cds_BDDb_im"};
        channelParameters[channel] = params;

        addSU3BreakingParameter("delta_A2t_cdcs_BDbD_re", "A2t_cscs_BDbD_re");
        addSU3BreakingParameter("delta_A2t_cdcs_BDbD_im", "A2t_cscs_BDbD_im");
        addAmplitudeParameter("dP3EW_ucs_BDbD_re", -ewp_limit, ewp_limit);
        addAmplitudeParameter("dP3EW_ucs_BDbD_im", -ewp_limit, ewp_limit);
        addAmplitudeParameter("R_A2_dcds_BDDb_re", -1., 1.1);
        addAmplitudeParameter("R_A2_dcds_BDDb_im", -1., 1.1);
        addSU3BreakingParameter("delta_G3t_cds_BDDb_re", "G3t_css_BDDb_re");
        addSU3BreakingParameter("delta_G3t_cds_BDDb_im", "G3t_css_BDDb_im");
    }
    else if (channel == "Bddspdsm")
    {
        // b → c(c̄d), spectator d
        vector<string> params = {
            "delta_A2t_cscd_BDbD_re", "delta_A2t_cscd_BDbD_im",
            "delta_G3t_csd_BDDb_re", "delta_G3t_csd_BDDb_im"};
        channelParameters[channel] = params;

        addSU3BreakingParameter("delta_A2t_cscd_BDbD_re", "A2t_cscs_BDbD_re");
        addSU3BreakingParameter("delta_A2t_cscd_BDbD_im", "A2t_cscs_BDbD_im");
        addSU3BreakingParameter("delta_G3t_csd_BDDb_re", "G3t_css_BDDb_re");
        addSU3BreakingParameter("delta_G3t_csd_BDDb_im", "G3t_css_BDDb_im");
    }
    else if (channel == "Bddspdm")
    {
        // b → c(c̄s), spectator d
        vector<string> params = {
            "delta_E1t_sccd_BDDb_re", "delta_E1t_sccd_BDDb_im",
            "delta_G1t_scd_BDDb_re", "delta_G1t_scd_BDDb_im"};
        channelParameters[channel] = params;

        addSU3BreakingParameter("delta_E1t_sccd_BDDb_re", "E1t_sccs_BDDb_re");
        addSU3BreakingParameter("delta_E1t_sccd_BDDb_im", "E1t_sccs_BDDb_im");
        addSU3BreakingParameter("delta_G1t_scd_BDDb_re", "G1t_scs_BDDb_re");
        addSU3BreakingParameter("delta_G1t_scd_BDDb_im", "G1t_scs_BDDb_im");
    }
    else if (channel == "Bddpdm")
    {
        // b → c(c̄d), spectator d
        vector<string> params = {
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
        vector<string> params = {
            "delta_A2t_cdcs_BDbD_re", "delta_A2t_cdcs_BDbD_im",
            "delta_A2t_cdcd_BDbD_re", "delta_A2t_cdcd_BDbD_im",
            "delta_dP3EW_ucd_BDbD_re", "delta_dP3EW_ucd_BDbD_im",
            "delta_A2_dcdd_BDDb_re", "delta_A2_dcdd_BDDb_im",
            "delta_G3t_cds_BDDb_re", "delta_G3t_cds_BDDb_im",
            "delta_G3t_cdd_BDDb_re", "delta_G3t_cdd_BDDb_im"};
        channelParameters[channel] = params;

        addSU3BreakingParameter("delta_A2t_cdcs_BDbD_re", "A2t_cscs_BDbD_re");
        addSU3BreakingParameter("delta_A2t_cdcs_BDbD_im", "A2t_cscs_BDbD_im");
        addSU3BreakingParameter("delta_A2t_cdcd_BDbD_re", "A2t_cdcs_BDbD_re");
        addSU3BreakingParameter("delta_A2t_cdcd_BDbD_im", "A2t_cdcs_BDbD_im");
        addSU3BreakingParameter("delta_dP3EW_ucd_BDbD_re", "dP3EW_ucs_BDbD_re");
        addSU3BreakingParameter("delta_dP3EW_ucd_BDbD_im", "dP3EW_ucs_BDbD_im");
        addSU3BreakingParameter("delta_A2_dcdd_BDDb_re", "A2_dcds_BDDb_re");
        addSU3BreakingParameter("delta_A2_dcdd_BDDb_im", "A2_dcds_BDDb_im");
        addSU3BreakingParameter("delta_G3t_cds_BDDb_re", "G3t_css_BDDb_re");
        addSU3BreakingParameter("delta_G3t_cds_BDDb_im", "G3t_css_BDDb_im");
        addSU3BreakingParameter("delta_G3t_cdd_BDDb_re", "G3t_cds_BDDb_re");
        addSU3BreakingParameter("delta_G3t_cdd_BDDb_im", "G3t_cds_BDDb_im");
    }
    else if (channel == "Bpdpd0b")
    {
        // b → c(c̄d), spectator u
        vector<string> params = {
            "delta_E1t_dccs_BDDb_re", "delta_E1t_dccs_BDDb_im",
            "delta_E1t_dccd_BDDb_re", "delta_E1t_dccd_BDDb_im",
            "dP1EW_dcu_BDDb_re", "dP1EW_dcu_BDDb_im",
            "R_A1_dcdd_BDDb_re", "R_A1_dcdd_BDDb_im",
            "delta_G1t_dcd_BDDb_re", "delta_G1t_dcd_BDDb_im",
            "delta_G1t_dcs_BDDb_re", "delta_G1t_dcs_BDDb_im"};
        channelParameters[channel] = params;

        addSU3BreakingParameter("delta_E1t_dccs_BDDb_re", "E1t_sccs_BDDb_re");
        addSU3BreakingParameter("delta_E1t_dccs_BDDb_im", "E1t_sccs_BDDb_im");
        addSU3BreakingParameter("delta_E1t_dccd_BDDb_re", "E1t_dccs_BDDb_re");
        addSU3BreakingParameter("delta_E1t_dccd_BDDb_im", "E1t_dccs_BDDb_im");
        addAmplitudeParameter("dP1EW_dcu_BDDb_re", -ewp_limit, ewp_limit);
        addAmplitudeParameter("dP1EW_dcu_BDDb_im", -ewp_limit, ewp_limit);
        addAmplitudeParameter("R_A1_dcdd_BDDb_re", -1., 1.1);
        addAmplitudeParameter("R_A1_dcdd_BDDb_im", -1., 1.1);
        addSU3BreakingParameter("delta_G1t_dcs_BDDb_re", "G1t_scs_BDDb_re");
        addSU3BreakingParameter("delta_G1t_dcs_BDDb_im", "G1t_scs_BDDb_im");
        addSU3BreakingParameter("delta_G1t_dcd_BDDb_re", "G1t_dcs_BDDb_re");
        addSU3BreakingParameter("delta_G1t_dcd_BDDb_im", "G1t_dcs_BDDb_im");
    }
    else if (channel == "Bpdspd0b")
    {
        // b → c(c̄s), spectator u
        vector<string> params = {
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
        cerr << "Error: Unrecognized channel \"" << channel << "\" in DefineParameters." << endl;
        throw runtime_error("Unrecognized channel: " + channel);
    }
}

map<string, double> parameterValues;
// function that given the channels and the string inside the map channelParameters adds every parameter name to a map <string, double> and returns said map.
map<string, double> goldenmodesB::DeclareParameters()
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
            cerr << "Channel " << channel << " not found in channelParameters." << endl;
        }
    }

    cout << "Declare Parameters called correctly" << endl;
    return parameterValues;
}

TComplex goldenmodesB::getPar(const string &baseName) const
{
    if (Debug)
    {
        cout << "Getting parameter: " << baseName << endl;
    }
    // Look for the real and imaginary parts in the parameterValues map
    auto it_real = parameterValues.find(baseName + "_re");
    auto it_imag = parameterValues.find(baseName + "_im");

    if (it_real != parameterValues.end() && it_imag != parameterValues.end())
    {
        // Parameter exists directly, return it
        if (Debug)
        {
            cout << "Found direct parameter: " << baseName << " = " << it_real->second << " + i*" << it_imag->second << endl;
        }
        return TComplex(it_real->second, it_imag->second);
    }

    // Check first if this is a ratio of amplitudes (has "R_" prefix)
    it_real = parameterValues.find("R_" + baseName + "_re");
    it_imag = parameterValues.find("R_" + baseName + "_im");

    if (it_real != parameterValues.end() && it_imag != parameterValues.end())
    {
        if (Debug)
        {
            cout << "Parameter " << baseName << " is a ratio with R = " << it_real->second << " + i*" << it_imag->second << endl;
        }
        // This is a ratio: result = reference * R
        TComplex R(it_real->second, it_imag->second);
        // Get the reference amplitude (handles chained ratios)
        string refName;
        if (baseName.find("JPSIP") != string::npos || baseName.find("PJPSI") != string::npos)
            refName = "E2t_ccsd_BJPSIP";
        else if (baseName.find("BJPSIV_0") != string::npos || baseName.find("BVJPSI_0") != string::npos)
            refName = "E2t_ccss_BJPSIV_0";
        else if (baseName.find("BJPSIV_paral") != string::npos || baseName.find("BVJPSI_paral") != string::npos)
        {
            refName = "E2t_ccss_BJPSIV_paral";
        }
        else if (baseName.find("BJPSIV_perp") != string::npos || baseName.find("BVJPSI_perp") != string::npos)
        {
            refName = "E2t_ccss_BJPSIV_perp";
        }
        else if (baseName.find("DDb") != string::npos || baseName.find("DbD") != string::npos)
        {
            refName = "E1t_sccs_BDDb";
        }
        else
        {
            cerr << "Error: Cannot determine base parameter for amplitude " << baseName << " so cannot define ratio parameter." << endl;
            exit(1);
        }
        if (Debug)
        {
            cout << "Parameter " << baseName << " is a ratio with reference " << refName << " and ratio R = " << R << endl;
            cout << "Recursively getting reference parameter " << refName << endl;
        }
        TComplex ref = getPar(refName); // Recursive call to get the reference amplitude
        if (Debug)
        {
            cout << "Reference amplitude " << refName << " = " << ref << endl;
            cout << "Calculating " << baseName << " = " << ref << " * " << R << endl;
        }
        return ref * R;
    }

    // Check if this is a derived amplitude (has a delta parameter)
    // Look for "<baseName>" in the deltaReferenceAmplitudes map
    auto deltaIt = deltaReferenceAmplitudes.find(baseName);
    if (deltaIt != deltaReferenceAmplitudes.end())
    {
        // This is a derived amplitude: result = reference * (1 + delta)
        // Get the delta value
        auto delta_real = parameterValues.find("delta_" + baseName + "_re");
        auto delta_imag = parameterValues.find("delta_" + baseName + "_im");

        if (delta_real != parameterValues.end() && delta_imag != parameterValues.end())
        {
            TComplex delta(delta_real->second, delta_imag->second);
            // Recursively get the reference amplitude (handles chained deltas)
            TComplex ref = getPar(deltaIt->second);
            if (Debug)
            {
                cout << "Parameter " << baseName << " is a derived amplitude with reference " << deltaIt->second << " and delta = " << delta << endl;
                cout << "Recursively getting reference parameter " << deltaIt->second << endl;
            }
            TComplex result = ref * (TComplex(1, 0) + delta);
            if (Debug)
            {
                cout << "Reference amplitude " << deltaIt->second << " = " << ref << endl;
                cout << "Calculating " << baseName << " = " << ref << " * (1 + " << delta << ") = " << result << endl;
            }
            return result;
        }
    }

    throw runtime_error("Error: Real or imaginary part for parameter " + baseName + " not found.");
}

// Setter function: sets the value for a given parameter in the map
void goldenmodesB::SetParameterValue(const string &paramName, double value)
{
    parameterValues[paramName] = value; // Insert or update the value for the given parameter name
}

//---------------------------------------------------------------------

double goldenmodesB::getParameterValue(const string &paramName) const
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
        cerr << "Error: TComplex " << paramName << " not found in parameterValues map." << endl;
        throw runtime_error("Parameter not found: " + paramName);
    }
}

//----------------------------------------------------------

// compute decay amplitude for each channel
void goldenmodesB::compute_decay_amplitudes(const string &channel)
{

    TComplex amp, ampc;
    TComplex amp_0, ampc_0;
    TComplex amp_paral, ampc_paral;
    TComplex amp_perp, ampc_perp;

    // Normalize the amplitudes by lam_bd_c for Bd mesons and by lam_bs_c for Bs mesons, so that we can factor out the mixing phase in the time-dependent CP asymmetries

    if (channel == "Bdjpsik0")
    {
        // Bd→J/ψ K⁰: b→c(c̄s), spectator d
        amp = lam_bs_c * getPar("E2t_ccsd_BJPSIP") - lam_bs_u * getPar("G2t_scd_BJPSIP");
        amp /= lam_bd_c; // Normalize by lam_bd_c
        ampc = lamst_bs_c * getPar("E2t_ccsd_BJPSIP") - lamst_bs_u * getPar("G2t_scd_BJPSIP");
        ampc /= lamst_bd_c; // Normalize by lamst_bd_c
        amplitude_map[channel] = make_pair(amp, ampc);
    }
    else if (channel == "Bdjpsip0")
    {
        // Bd→J/ψ π⁰: b→c(c̄d), spectator d
        amp = (lam_bd_c * (getPar("E2t_ccdd_BJPSIP") + getPar("dP4EW_ucd_BPJPSI")) -
               lam_bd_u * (getPar("EA2_ddcd_BPJPSI") + getPar("dP4EW_ucd_BPJPSI") + getPar("G2t_dcd_BJPSIP"))) /
              sqrt(2.);
        amp /= lam_bd_c; // Normalize by lam_bd_c
        ampc = (lamst_bd_c * (getPar("E2t_ccdd_BJPSIP") + getPar("dP4EW_ucd_BPJPSI")) -
                lamst_bd_u * (getPar("EA2_ddcd_BPJPSI") + getPar("dP4EW_ucd_BPJPSI") + getPar("G2t_dcd_BJPSIP"))) /
               sqrt(2.);
        ampc /= lamst_bd_c; // Normalize by lamst_bd_c
        amplitude_map[channel] = make_pair(amp, ampc);
    }
    else if (channel == "Bdjpsieta8")
    {
        // Bd→J/ψ eta_8: b→c(c̄d), spectator d
        amp = ((lam_bd_c * (getPar("E2t_ccdd_BJPSIP") + getPar("dP4EW_ucd_BPJPSI") + 2. * getPar("EA2t_ccdd_BJPSIP") - 2. * getPar("EA2t_ccsd_BJPSIP"))) +
               lam_bd_u * (getPar("EA2_ddcd_BPJPSI") + getPar("dP4EW_ucd_BPJPSI") - getPar("G2t_dcd_BJPSIP") - 2. * getPar("G4t_cdd_BJPSIP") + 2. * getPar("G4t_csd_BJPSIP"))) /
              sqrt(6.);
        amp /= lam_bd_c; // Normalize by lam_bd_c
        ampc = ((lamst_bd_c * (getPar("E2t_ccdd_BJPSIP") + getPar("dP4EW_ucd_BPJPSI") + 2. * getPar("EA2t_ccdd_BJPSIP") - 2. * getPar("EA2t_ccsd_BJPSIP"))) +
                lamst_bd_u * (getPar("EA2_ddcd_BPJPSI") + getPar("dP4EW_ucd_BPJPSI") - getPar("G2t_dcd_BJPSIP") - 2. * getPar("G4t_cdd_BJPSIP") + 2. * getPar("G4t_csd_BJPSIP"))) /
               sqrt(6.);
        ampc = amp / lamst_bd_c; // Normalize by lamst_bd_c
        amplitude_map[channel] = make_pair(amp, ampc);
    }
    else if (channel == "Bdjpsieta1")
    {
        // Bd→J/ψ eta_1: b→c(c̄d), spectator d
        amp = ((lam_bd_c * (getPar("E2t_ccdd_BJPSIP") + getPar("dP4EW_ucd_BPJPSI") + 2. * getPar("EA2t_ccdd_BJPSIP") + getPar("EA2t_ccsd_BJPSIP"))) +
               lam_bd_u * (getPar("EA2_ddcd_BPJPSI") + getPar("dP4EW_ucd_BPJPSI") - getPar("G2t_dcd_BJPSIP") - 2. * getPar("G4t_cdd_BJPSIP") - getPar("G4t_csd_BJPSIP"))) /
              sqrt(3.);
        amp /= lam_bd_c; // Normalize by lam_bd_c
        ampc = ((lamst_bd_c * (getPar("E2t_ccdd_BJPSIP") + getPar("dP4EW_ucd_BPJPSI") + 2. * getPar("EA2t_ccdd_BJPSIP") + getPar("EA2t_ccsd_BJPSIP"))) +
                lamst_bd_u * (getPar("EA2_ddcd_BPJPSI") + getPar("dP4EW_ucd_BPJPSI") - getPar("G2t_dcd_BJPSIP") - 2. * getPar("G4t_cdd_BJPSIP") - getPar("G4t_csd_BJPSIP"))) /
               sqrt(3.);
        ampc /= lamst_bd_c; // Normalize by lamst_bd_c
        amplitude_map[channel] = make_pair(amp, ampc);
    }
    else if (channel == "Bpjpsikp")
    {
        // B⁺→J/ψ K⁺: b→c(c̄s), spectator u
        amp = lam_bs_c * (getPar("E2t_ccsd_BJPSIP") + getPar("dP2EW_scu_BPJPSI")) +
              lam_bs_u * (getPar("EA1_sdcd_BPJPSI") + getPar("dP2EW_scu_BPJPSI") - getPar("G2t_scd_BJPSIP"));
        ampc = lamst_bs_c * (getPar("E2t_ccsd_BJPSIP") + getPar("dP2EW_scu_BPJPSI")) +
               lamst_bs_u * (getPar("EA1_sdcd_BPJPSI") + getPar("dP2EW_scu_BPJPSI") - getPar("G2t_scd_BJPSIP"));
        amplitude_map[channel] = make_pair(amp, ampc);
    }
    else if (channel == "Bpjpsipp")
    {
        // B⁺→J/ψ π⁺: b→c(c̄d), spectator u
        amp = lam_bd_c * (getPar("E2t_ccdd_BJPSIP") + getPar("dP2EW_dcu_BPJPSI")) +
              lam_bd_u * (getPar("EA1_ddcd_BPJPSI") + getPar("dP2EW_dcu_BPJPSI") - getPar("G2t_dcd_BJPSIP"));
        ampc = lamst_bd_c * (getPar("E2t_ccdd_BJPSIP") + getPar("dP2EW_dcu_BPJPSI")) +
               lamst_bd_u * (getPar("EA1_ddcd_BPJPSI") + getPar("dP2EW_dcu_BPJPSI") - getPar("G2t_dcd_BJPSIP"));
        amplitude_map[channel] = make_pair(amp, ampc);
    }
    else if (channel == "Bsjpsip0")
    {
        // Bs→J/ψ π⁰: b→c(c̄s), spectator s
        amp = -(lam_bs_c * getPar("dP4EW_ucs_BPJPSI") + lam_bs_u * (getPar("EA2_ddcs_BPJPSI") + getPar("dP4EW_ucs_BPJPSI"))) / sqrt(2.);
        amp /= lam_bs_c; // Normalize by lam_bs_c
        ampc = -(lamst_bs_c * getPar("dP4EW_ucs_BPJPSI") + lamst_bs_u * (getPar("EA2_ddcs_BPJPSI") + getPar("dP4EW_ucs_BPJPSI"))) / sqrt(2.);
        ampc /= lamst_bs_c; // Normalize by lamst_bs_c
        amplitude_map[channel] = make_pair(amp, ampc);
    }
    else if (channel == "Bsjpsik0b")
    {
        // Bs→J/ψ \bar{K}⁰: b→c(c̄d), spectator s
        amp = lam_bd_c * getPar("E2t_ccds_BJPSIP") - lam_bd_u * getPar("G2t_dcs_BJPSIP");
        amp /= lam_bs_c; // Normalize by lam_bs_c
        ampc = lamst_bd_c * getPar("E2t_ccds_BJPSIP") - lamst_bd_u * getPar("G2t_dcs_BJPSIP");
        ampc /= lamst_bs_c; // Normalize by lamst_bs_c
        amplitude_map[channel] = make_pair(amp, ampc);
    }
    else if (channel == "Bsjpsieta8")
    {
        // Bs→J/ψ eta_8: b→c(c̄s), spectator s
        amp = ((lam_bs_c * (-2. * getPar("E2t_ccss_BJPSIP") + getPar("dP4EW_ucs_BPJPSI") + 2. * getPar("EA2t_ccds_BJPSIP") - 2. * getPar("EA2t_ccss_BJPSIP"))) +
               lam_bs_u * (getPar("EA2_ddcs_BPJPSI") + getPar("dP4EW_ucs_BPJPSI") + 2. * getPar("G2t_scs_BJPSIP") - 2. * getPar("G4t_cds_BJPSIP") + 2. * getPar("G4t_css_BJPSIP"))) /
              sqrt(6.);
        amp /= lam_bs_c; // Normalize by lam_bs_c
        ampc = ((lamst_bs_c * (-2. * getPar("E2t_ccss_BJPSIP") + getPar("dP4EW_ucs_BPJPSI") + 2. * getPar("EA2t_ccds_BJPSIP") - 2. * getPar("EA2t_ccss_BJPSIP"))) +
                lamst_bs_u * (getPar("EA2_ddcs_BPJPSI") + getPar("dP4EW_ucs_BPJPSI") + 2. * getPar("G2t_scs_BJPSIP") - 2. * getPar("G4t_cds_BJPSIP") + 2. * getPar("G4t_css_BJPSIP"))) /
               sqrt(6.);
        ampc /= lamst_bs_c; // Normalize by lamst_bs_c
        amplitude_map[channel] = make_pair(amp, ampc);
    }
    else if (channel == "Bsjpsieta1")
    {
        // Bs→J/ψ eta_1: b→c(c̄s), spectator s
        amp = ((lam_bs_c * (getPar("E2t_ccss_BJPSIP") + getPar("dP4EW_ucs_BPJPSI") + 2. * getPar("EA2t_ccds_BJPSIP") + getPar("EA2t_ccss_BJPSIP"))) +
               lam_bs_u * (getPar("EA2_ddcs_BPJPSI") + getPar("dP4EW_ucs_BPJPSI") - getPar("G2t_scs_BJPSIP") - 2. * getPar("G4t_cds_BJPSIP") - getPar("G4t_css_BJPSIP"))) /
              sqrt(3.);
        amp /= lam_bs_c; // Normalize by lam_bs_c
        ampc = ((lamst_bs_c * (getPar("E2t_ccss_BJPSIP") + getPar("dP4EW_ucs_BPJPSI") + 2. * getPar("EA2t_ccds_BJPSIP") + getPar("EA2t_ccss_BJPSIP"))) +
                lamst_bs_u * (getPar("EA2_ddcs_BPJPSI") + getPar("dP4EW_ucs_BPJPSI") - getPar("G2t_scs_BJPSIP") - 2. * getPar("G4t_cds_BJPSIP") - getPar("G4t_css_BJPSIP"))) /
               sqrt(3.);
        ampc /= lamst_bs_c; // Normalize by lamst_bs_c
        amplitude_map[channel] = make_pair(amp, ampc);
    }
    else if (channel == "Bsjpsiphi")
    {
        // Bs→J/ψ φ: b→c(c̄s), spectator s
        amp_0 = -lam_bs_c * (getPar("E2t_ccss_BJPSIV_0") + getPar("EA2t_ccss_BJPSIV_0")) +
                lam_bs_u * (getPar("G2t_scs_BJPSIV_0") + getPar("G4t_css_BJPSIV_0"));
        amp_0 /= lam_bs_c; // Normalize by lam_bs_c
        ampc_0 = -lamst_bs_c * (getPar("E2t_ccss_BJPSIV_0") + getPar("EA2t_ccss_BJPSIV_0")) +
                 lamst_bs_u * (getPar("G2t_scs_BJPSIV_0") + getPar("G4t_css_BJPSIV_0"));
        ampc_0 /= lamst_bs_c; // Normalize by lamst_bs_c
        amp_paral = -lam_bs_c * (getPar("E2t_ccss_BJPSIV_paral") + getPar("EA2t_ccss_BJPSIV_paral")) +
                    lam_bs_u * (getPar("G2t_scs_BJPSIV_paral") + getPar("G4t_css_BJPSIV_paral"));
        amp_paral /= lam_bs_c; // Normalize by lam_bs_c
        ampc_paral = -lamst_bs_c * (getPar("E2t_ccss_BJPSIV_paral") + getPar("EA2t_ccss_BJPSIV_paral")) +
                     lamst_bs_u * (getPar("G2t_scs_BJPSIV_paral") + getPar("G4t_css_BJPSIV_paral"));
        ampc_paral /= lamst_bs_c; // Normalize by lamst_bs_c
        amp_perp = -lam_bs_c * (getPar("E2t_ccss_BJPSIV_perp") + getPar("EA2t_ccss_BJPSIV_perp")) +
                   lam_bs_u * (getPar("G2t_scs_BJPSIV_perp") + getPar("G4t_css_BJPSIV_perp"));
        amp_perp /= lam_bs_c; // Normalize by lam_bs_c
        ampc_perp = -lamst_bs_c * (getPar("E2t_ccss_BJPSIV_perp") + getPar("EA2t_ccss_BJPSIV_perp")) +
                    lamst_bs_u * (getPar("G2t_scs_BJPSIV_perp") + getPar("G4t_css_BJPSIV_perp"));
        ampc_perp /= lamst_bs_c; // Normalize by lamst_bs_c
        amplitude_map[channel + "_0"] = make_pair(amp_0, ampc_0);
        amplitude_map[channel + "_paral"] = make_pair(amp_paral, ampc_paral);
        amplitude_map[channel + "_perp"] = make_pair(amp_perp, ampc_perp);
    }
    else if (channel == "Bsjpsiom")
    {
        // Bs→J/ψ ω: b→c(c̄s), spectator s
        amp_0 = (lam_bs_c * (2. * getPar("EA2t_ccds_BJPSIV_0") + getPar("dP4EW_ucs_BVJPSI_0")) +
                 lam_bs_u * (getPar("EA2_ddcs_BVJPSI_0") + getPar("dP4EW_ucs_BVJPSI_0") - 2. * getPar("G4t_cds_BJPSIV_0"))) /
                sqrt(2.);
        amp_0 /= lam_bs_c; // Normalize by lam_bs_c
        ampc_0 = (lamst_bs_c * (2. * getPar("EA2t_ccds_BJPSIV_0") + getPar("dP4EW_ucs_BVJPSI_0")) +
                  lamst_bs_u * (getPar("EA2_ddcs_BVJPSI_0") + getPar("dP4EW_ucs_BVJPSI_0") - 2. * getPar("G4t_cds_BJPSIV_0"))) /
                 sqrt(2.);
        ampc_0 /= lamst_bs_c; // Normalize by lamst_bs_c
        amp_paral = (lam_bs_c * (2. * getPar("EA2t_ccds_BJPSIV_paral") + getPar("dP4EW_ucs_BVJPSI_paral")) +
                     lam_bs_u * (getPar("EA2_ddcs_BVJPSI_paral") + getPar("dP4EW_ucs_BVJPSI_paral") - 2. * getPar("G4t_cds_BJPSIV_paral"))) /
                    sqrt(2.);
        amp_paral /= lam_bs_c; // Normalize by lam_bs_c
        ampc_paral = (lamst_bs_c * (2. * getPar("EA2t_ccds_BJPSIV_paral") + getPar("dP4EW_ucs_BVJPSI_paral")) +
                      lamst_bs_u * (getPar("EA2_ddcs_BVJPSI_paral") + getPar("dP4EW_ucs_BVJPSI_paral") - 2. * getPar("G4t_cds_BJPSIV_paral"))) /
                     sqrt(2.);
        ampc_paral /= lamst_bs_c; // Normalize by lamst_bs_c
        amp_perp = (lam_bs_c * (2. * getPar("EA2t_ccds_BJPSIV_perp") + getPar("dP4EW_ucs_BVJPSI_perp")) +
                    lam_bs_u * (getPar("EA2_ddcs_BVJPSI_perp") + getPar("dP4EW_ucs_BVJPSI_perp") - 2. * getPar("G4t_cds_BJPSIV_perp"))) /
                   sqrt(2.);
        amp_perp /= lam_bs_c; // Normalize by lam_bs_c
        ampc_perp = (lamst_bs_c * (2. * getPar("EA2t_ccds_BJPSIV_perp") + getPar("dP4EW_ucs_BVJPSI_perp")) +
                     lamst_bs_u * (getPar("EA2_ddcs_BVJPSI_perp") + getPar("dP4EW_ucs_BVJPSI_perp") - 2. * getPar("G4t_cds_BJPSIV_perp"))) /
                    sqrt(2.);
        ampc_perp /= lamst_bs_c; // Normalize by lamst_bs_c
        amplitude_map[channel + "_0"] = make_pair(amp_0, ampc_0);
        amplitude_map[channel + "_paral"] = make_pair(amp_paral, ampc_paral);
        amplitude_map[channel + "_perp"] = make_pair(amp_perp, ampc_perp);
    }
    else if (channel == "Bsjpsikbst0")
    {
        // Bs→J/ψ \bar{K}*: b→c(c̄d), spectator s
        amp_0 = -lam_bd_c * getPar("E2t_ccds_BJPSIV_0") - lam_bd_u * getPar("G2t_dcs_BJPSIV_0");
        amp_0 /= lam_bs_c; // Normalize by lam_bs_c
        ampc_0 = -lamst_bd_c * getPar("E2t_ccds_BJPSIV_0") - lamst_bd_u * getPar("G2t_dcs_BJPSIV_0");
        ampc_0 /= lamst_bs_c; // Normalize by lamst_bs_c
        amp_paral = -lam_bd_c * getPar("E2t_ccds_BJPSIV_paral") - lam_bd_u * getPar("G2t_dcs_BJPSIV_paral");
        amp_paral /= lam_bs_c; // Normalize by lam_bs_c
        ampc_paral = -lamst_bd_c * getPar("E2t_ccds_BJPSIV_paral") - lamst_bd_u * getPar("G2t_dcs_BJPSIV_paral");
        ampc_paral /= lamst_bs_c; // Normalize by lamst_bs_c
        amp_perp = -lam_bd_c * getPar("E2t_ccds_BJPSIV_perp") - lam_bd_u * getPar("G2t_dcs_BJPSIV_perp");
        amp_perp /= lam_bs_c; // Normalize by lam_bs_c
        ampc_perp = -lamst_bd_c * getPar("E2t_ccds_BJPSIV_perp") - lamst_bd_u * getPar("G2t_dcs_BJPSIV_perp");
        ampc_perp /= lamst_bs_c; // Normalize by lamst_bs_c
        amplitude_map[channel + "_0"] = make_pair(amp_0, ampc_0);
        amplitude_map[channel + "_paral"] = make_pair(amp_paral, ampc_paral);
        amplitude_map[channel + "_perp"] = make_pair(amp_perp, ampc_perp);
    }
    else if (channel == "Bsjpsirho0")
    {
        // Bs→J/ψ \rho⁰: b→c(c̄s), spectator s
        amp_0 = -(lam_bs_c * getPar("dP4EW_ucs_BVJPSI_0") + lam_bs_u * (getPar("EA2_ddcs_BVJPSI_0") + getPar("dP4EW_ucs_BVJPSI_0"))) / sqrt(2.);
        amp_0 /= lam_bs_c; // Normalize by lam_bs_c
        ampc_0 = -(lamst_bs_c * getPar("dP4EW_ucs_BVJPSI_0") + lamst_bs_u * (getPar("EA2_ddcs_BVJPSI_0") + getPar("dP4EW_ucs_BVJPSI_0"))) / sqrt(2.);
        ampc_0 /= lamst_bs_c; // Normalize by lamst_bs_c
        amp_paral = -(lam_bs_c * getPar("dP4EW_ucs_BVJPSI_paral") + lam_bs_u * (getPar("EA2_ddcs_BVJPSI_paral") + getPar("dP4EW_ucs_BVJPSI_paral"))) / sqrt(2.);
        amp_paral /= lam_bs_c; // Normalize by lam_bs_c
        ampc_paral = -(lamst_bs_c * getPar("dP4EW_ucs_BVJPSI_paral") + lamst_bs_u * (getPar("EA2_ddcs_BVJPSI_paral") + getPar("dP4EW_ucs_BVJPSI_paral"))) / sqrt(2.);
        ampc_paral /= lamst_bs_c; // Normalize by lamst_bs_c
        amp_perp = -(lam_bs_c * getPar("dP4EW_ucs_BVJPSI_perp") + lam_bs_u * (getPar("EA2_ddcs_BVJPSI_perp") + getPar("dP4EW_ucs_BVJPSI_perp"))) / sqrt(2.);
        amp_perp /= lam_bs_c; // Normalize by lam_bs_c
        ampc_perp = -(lamst_bs_c * getPar("dP4EW_ucs_BVJPSI_perp") + lamst_bs_u * (getPar("EA2_ddcs_BVJPSI_perp") + getPar("dP4EW_ucs_BVJPSI_perp"))) / sqrt(2.);
        ampc_perp /= lamst_bs_c; // Normalize by lamst_bs_c
        amplitude_map[channel + "_0"] = make_pair(amp_0, ampc_0);
        amplitude_map[channel + "_paral"] = make_pair(amp_paral, ampc_paral);
        amplitude_map[channel + "_perp"] = make_pair(amp_perp, ampc_perp);
    }
    else if (channel == "Bdjpsiom")
    {
        // Bd→J/ψ ω: b→c(c̄d), spectator d
        amp_0 = (lam_bd_c * (getPar("E2t_ccdd_BJPSIV_0") + 2. * getPar("EA2t_ccdd_BJPSIV_0") + getPar("dP4EW_ucd_BVJPSI_0")) +
                 lam_bd_u * (getPar("EA2_ddcd_BVJPSI_0") + getPar("dP4EW_ucd_BVJPSI_0") - getPar("G2t_dcd_BJPSIV_0") - 2. * getPar("G4t_cdd_BJPSIV_0"))) /
                sqrt(2.);
        amp_0 /= lam_bd_c; // Normalize by lam_bd_c
        ampc_0 = (lamst_bd_c * (getPar("E2t_ccdd_BJPSIV_0") + 2. * getPar("EA2t_ccdd_BJPSIV_0") + getPar("dP4EW_ucd_BVJPSI_0")) +
                  lamst_bd_u * (getPar("EA2_ddcd_BVJPSI_0") + getPar("dP4EW_ucd_BVJPSI_0") - getPar("G2t_dcd_BJPSIV_0") - 2. * getPar("G4t_cdd_BJPSIV_0"))) /
                 sqrt(2.);
        ampc_0 /= lamst_bd_c; // Normalize by lamst_bd_c
        amp_paral = (lam_bd_c * (getPar("E2t_ccdd_BJPSIV_paral") + 2. * getPar("EA2t_ccdd_BJPSIV_paral") + getPar("dP4EW_ucd_BVJPSI_paral")) +
                     lam_bd_u * (getPar("EA2_ddcd_BVJPSI_paral") + getPar("dP4EW_ucd_BVJPSI_paral") - getPar("G2t_dcd_BJPSIV_paral") - 2. * getPar("G4t_cdd_BJPSIV_paral"))) /
                    sqrt(2.);
        amp_paral /= lam_bd_c; // Normalize by lam_bd_c
        ampc_paral = (lamst_bd_c * (getPar("E2t_ccdd_BJPSIV_paral") + 2. * getPar("EA2t_ccdd_BJPSIV_paral") + getPar("dP4EW_ucd_BVJPSI_paral")) +
                      lamst_bd_u * (getPar("EA2_ddcd_BVJPSI_paral") + getPar("dP4EW_ucd_BVJPSI_paral") - getPar("G2t_dcd_BJPSIV_paral") - 2. * getPar("G4t_cdd_BJPSIV_paral"))) /
                     sqrt(2.);
        ampc_paral /= lamst_bd_c; // Normalize by lamst_bd_c
        amp_perp = (lam_bd_c * (getPar("E2t_ccdd_BJPSIV_perp") + 2. * getPar("EA2t_ccdd_BJPSIV_perp") + getPar("dP4EW_ucd_BVJPSI_perp")) +
                    lam_bd_u * (getPar("EA2_ddcd_BVJPSI_perp") + getPar("dP4EW_ucd_BVJPSI_perp") - getPar("G2t_dcd_BJPSIV_perp") - 2. * getPar("G4t_cdd_BJPSIV_perp"))) /
                   sqrt(2.);
        amp_perp /= lam_bd_c; // Normalize by lam_bd_c
        ampc_perp = (lamst_bd_c * (getPar("E2t_ccdd_BJPSIV_perp") + 2. * getPar("EA2t_ccdd_BJPSIV_perp") + getPar("dP4EW_ucd_BVJPSI_perp")) +
                     lamst_bd_u * (getPar("EA2_ddcd_BVJPSI_perp") + getPar("dP4EW_ucd_BVJPSI_perp") - getPar("G2t_dcd_BJPSIV_perp") - 2. * getPar("G4t_cdd_BJPSIV_perp"))) /
                    sqrt(2.);
        ampc_perp /= lamst_bd_c; // Normalize by lamst_bd_c
        amplitude_map[channel + "_0"] = make_pair(amp_0, ampc_0);
        amplitude_map[channel + "_paral"] = make_pair(amp_paral, ampc_paral);
        amplitude_map[channel + "_perp"] = make_pair(amp_perp, ampc_perp);
    }
    else if (channel == "Bdjpsikst0")
    {
        // Bd→J/ψ K*: b→c(c̄s), spectator d
        amp_0 = lam_bs_c * getPar("E2t_ccsd_BJPSIV_0") - lam_bs_u * getPar("G2t_scd_BJPSIV_0");
        amp_0 /= lam_bd_c; // Normalize by lam_bd_c
        ampc_0 = lamst_bs_c * getPar("E2t_ccsd_BJPSIV_0") - lamst_bs_u * getPar("G2t_scd_BJPSIV_0");
        ampc_0 /= lamst_bd_c; // Normalize by lamst_bd_c
        amp_paral = lam_bs_c * getPar("E2t_ccsd_BJPSIV_paral") - lam_bs_u * getPar("G2t_scd_BJPSIV_paral");
        amp_paral /= lam_bd_c; // Normalize by lam_bd_c
        ampc_paral = lamst_bs_c * getPar("E2t_ccsd_BJPSIV_paral") - lamst_bs_u * getPar("G2t_scd_BJPSIV_paral");
        ampc_paral /= lamst_bd_c; // Normalize by lamst_bd_c
        amp_perp = lam_bs_c * getPar("E2t_ccsd_BJPSIV_perp") - lam_bs_u * getPar("G2t_scd_BJPSIV_perp");
        amp_perp /= lam_bd_c; // Normalize by lam_bd_c
        ampc_perp = lamst_bs_c * getPar("E2t_ccsd_BJPSIV_perp") - lamst_bs_u * getPar("G2t_scd_BJPSIV_perp");
        ampc_perp /= lamst_bd_c; // Normalize by lamst_bd_c
        amplitude_map[channel + "_0"] = make_pair(amp_0, ampc_0);
        amplitude_map[channel + "_paral"] = make_pair(amp_paral, ampc_paral);
        amplitude_map[channel + "_perp"] = make_pair(amp_perp, ampc_perp);
    }
    else if (channel == "Bdjpsirho0")
    {
        // Bd→J/ψ ρ: b→c(c̄d), spectator d
        amp_0 = (lam_bd_c * (getPar("E2t_ccdd_BJPSIV_0") + getPar("dP4EW_ucd_BVJPSI_0")) -
                 lam_bd_u * (getPar("EA2_ddcd_BVJPSI_0") + getPar("dP4EW_ucd_BVJPSI_0") + getPar("G2t_dcd_BJPSIV_0"))) /
                sqrt(2.);
        amp_0 /= lam_bd_c; // Normalize by lam_bd_c
        ampc_0 = (lamst_bd_c * (getPar("E2t_ccdd_BJPSIV_0") + getPar("dP4EW_ucd_BVJPSI_0")) -
                  lamst_bd_u * (getPar("EA2_ddcd_BVJPSI_0") + getPar("dP4EW_ucd_BVJPSI_0") + getPar("G2t_dcd_BJPSIV_0"))) /
                 sqrt(2.);
        ampc_0 /= lamst_bd_c; // Normalize by lamst_bd_c
        amp_paral = (lam_bd_c * (getPar("E2t_ccdd_BJPSIV_paral") + getPar("dP4EW_ucd_BVJPSI_paral")) -
                     lam_bd_u * (getPar("EA2_ddcd_BVJPSI_paral") + getPar("dP4EW_ucd_BVJPSI_paral") + getPar("G2t_dcd_BJPSIV_paral"))) /
                    sqrt(2.);
        amp_paral /= lam_bd_c; // Normalize by lam_bd_c
        ampc_paral = (lamst_bd_c * (getPar("E2t_ccdd_BJPSIV_paral") + getPar("dP4EW_ucd_BVJPSI_paral")) -
                      lamst_bd_u * (getPar("EA2_ddcd_BVJPSI_paral") + getPar("dP4EW_ucd_BVJPSI_paral") + getPar("G2t_dcd_BJPSIV_paral"))) /
                     sqrt(2.);
        ampc_paral /= lamst_bd_c; // Normalize by lamst_bd_c
        amp_perp = (lam_bd_c * (getPar("E2t_ccdd_BJPSIV_perp") + getPar("dP4EW_ucd_BVJPSI_perp")) -
                    lam_bd_u * (getPar("EA2_ddcd_BVJPSI_perp") + getPar("dP4EW_ucd_BVJPSI_perp") + getPar("G2t_dcd_BJPSIV_perp"))) /
                   sqrt(2.);
        amp_perp /= lam_bd_c; // Normalize by lam_bd_c
        ampc_perp = (lamst_bd_c * (getPar("E2t_ccdd_BJPSIV_perp") + getPar("dP4EW_ucd_BVJPSI_perp")) -
                     lamst_bd_u * (getPar("EA2_ddcd_BVJPSI_perp") + getPar("dP4EW_ucd_BVJPSI_perp") + getPar("G2t_dcd_BJPSIV_perp"))) /
                    sqrt(2.);
        ampc_perp /= lamst_bd_c; // Normalize by lamst_bd_c
        amplitude_map[channel + "_0"] = make_pair(amp_0, ampc_0);
        amplitude_map[channel + "_paral"] = make_pair(amp_paral, ampc_paral);
        amplitude_map[channel + "_perp"] = make_pair(amp_perp, ampc_perp);
    }
    else if (channel == "Bdjpsiphi")
    {
        // Bd→J/ψ \phi: b→c(c̄d), spectator d
        amp_0 = -lam_bd_c * getPar("EA2t_ccsd_BJPSIV_0") - lam_bd_u * getPar("G4t_csd_BJPSIV_0");
        ampc_0 = -lamst_bd_c * getPar("EA2t_ccsd_BJPSIV_0") - lamst_bd_u * getPar("G4t_csd_BJPSIV_0");
        amp_paral = -lam_bd_c * getPar("EA2t_ccsd_BJPSIV_paral") - lam_bd_u * getPar("G4t_csd_BJPSIV_paral");
        ampc_paral = -lamst_bd_c * getPar("EA2t_ccsd_BJPSIV_paral") - lamst_bd_u * getPar("G4t_csd_BJPSIV_paral");
        amp_perp = -lam_bd_c * getPar("EA2t_ccsd_BJPSIV_perp") - lam_bd_u * getPar("G4t_csd_BJPSIV_perp");
        ampc_perp = -lamst_bd_c * getPar("EA2t_ccsd_BJPSIV_perp") - lamst_bd_u * getPar("G4t_csd_BJPSIV_perp");
        amp_0 /= lam_bd_c;        // Normalize by lam_bd_c
        ampc_0 /= lamst_bd_c;     // Normalize by lamst_bd_c
        amp_paral /= lam_bd_c;    // Normalize by lam_bd_c
        ampc_paral /= lamst_bd_c; // Normalize by lamst_bd_c
        amp_perp /= lam_bd_c;     // Normalize by lam_bd_c
        ampc_perp /= lamst_bd_c;  // Normalize by lamst_bd_c
        amplitude_map[channel + "_0"] = make_pair(amp_0, ampc_0);
        amplitude_map[channel + "_paral"] = make_pair(amp_paral, ampc_paral);
        amplitude_map[channel + "_perp"] = make_pair(amp_perp, ampc_perp);
    }
    else if (channel == "Bpjpsirhop")
    {
        // B⁺→J/ψ ρ⁺: b→c(c̄d), spectator u
        amp_0 = lam_bd_c * (getPar("E2t_ccdd_BJPSIV_0") + getPar("dP2EW_dcu_BJPSIV_0")) +
                lam_bd_u * (getPar("EA1_ddcd_BVJPSI_0") + getPar("dP2EW_dcu_BJPSIV_0") - getPar("G2t_dcd_BJPSIV_0"));
        ampc_0 = lamst_bd_c * (getPar("E2t_ccdd_BJPSIV_0") + getPar("dP2EW_dcu_BJPSIV_0")) +
                 lamst_bd_u * (getPar("EA1_ddcd_BVJPSI_0") + getPar("dP2EW_dcu_BJPSIV_0") - getPar("G2t_dcd_BJPSIV_0"));
        amp_paral = lam_bd_c * (getPar("E2t_ccdd_BJPSIV_paral") + getPar("dP2EW_dcu_BJPSIV_paral")) +
                    lam_bd_u * (getPar("EA1_ddcd_BVJPSI_paral") + getPar("dP2EW_dcu_BJPSIV_paral") - getPar("G2t_dcd_BJPSIV_paral"));
        ampc_paral = lamst_bd_c * (getPar("E2t_ccdd_BJPSIV_paral") + getPar("dP2EW_dcu_BJPSIV_paral")) +
                     lamst_bd_u * (getPar("EA1_ddcd_BVJPSI_paral") + getPar("dP2EW_dcu_BJPSIV_paral") - getPar("G2t_dcd_BJPSIV_paral"));
        amp_perp = lam_bd_c * (getPar("E2t_ccdd_BJPSIV_perp") + getPar("dP2EW_dcu_BJPSIV_perp")) +
                   lam_bd_u * (getPar("EA1_ddcd_BVJPSI_perp") + getPar("dP2EW_dcu_BJPSIV_perp") - getPar("G2t_dcd_BJPSIV_perp"));
        ampc_perp = lamst_bd_c * (getPar("E2t_ccdd_BJPSIV_perp") + getPar("dP2EW_dcu_BJPSIV_perp")) +
                    lamst_bd_u * (getPar("EA1_ddcd_BVJPSI_perp") + getPar("dP2EW_dcu_BJPSIV_perp") - getPar("G2t_dcd_BJPSIV_perp"));
        amplitude_map[channel + "_0"] = make_pair(amp_0, ampc_0);
        amplitude_map[channel + "_paral"] = make_pair(amp_paral, ampc_paral);
        amplitude_map[channel + "_perp"] = make_pair(amp_perp, ampc_perp);
    }
    else if (channel == "Bpjpsikstp")
    {
        // B⁺→J/ψ K*⁺: b→c(c̄s), spectator u
        amp_0 = lam_bs_c * (getPar("E2t_ccsd_BJPSIV_0") + getPar("dP2EW_scu_BJPSIV_0")) +
                lam_bs_u * (getPar("EA1_sdcd_BVJPSI_0") + getPar("dP2EW_scu_BJPSIV_0") - getPar("G2t_scd_BJPSIV_0"));
        ampc_0 = lamst_bs_c * (getPar("E2t_ccsd_BJPSIV_0") + getPar("dP2EW_scu_BJPSIV_0")) +
                 lamst_bs_u * (getPar("EA1_sdcd_BVJPSI_0") + getPar("dP2EW_scu_BJPSIV_0") - getPar("G2t_scd_BJPSIV_0"));
        amp_paral = lam_bs_c * (getPar("E2t_ccsd_BJPSIV_paral") + getPar("dP2EW_scu_BJPSIV_paral")) +
                    lam_bs_u * (getPar("EA1_sdcd_BVJPSI_paral") + getPar("dP2EW_scu_BJPSIV_paral") - getPar("G2t_scd_BJPSIV_paral"));
        ampc_paral = lamst_bs_c * (getPar("E2t_ccsd_BJPSIV_paral") + getPar("dP2EW_scu_BJPSIV_paral")) +
                     lamst_bs_u * (getPar("EA1_sdcd_BVJPSI_paral") + getPar("dP2EW_scu_BJPSIV_paral") - getPar("G2t_scd_BJPSIV_paral"));
        amp_perp = lam_bs_c * (getPar("E2t_ccsd_BJPSIV_perp") + getPar("dP2EW_scu_BJPSIV_perp")) +
                   lam_bs_u * (getPar("EA1_sdcd_BVJPSI_perp") + getPar("dP2EW_scu_BJPSIV_perp") - getPar("G2t_scd_BJPSIV_perp"));
        ampc_perp = lamst_bs_c * (getPar("E2t_ccsd_BJPSIV_perp") + getPar("dP2EW_scu_BJPSIV_perp")) +
                    lamst_bs_u * (getPar("EA1_sdcd_BVJPSI_perp") + getPar("dP2EW_scu_BJPSIV_perp") - getPar("G2t_scd_BJPSIV_perp"));
        amplitude_map[channel + "_0"] = make_pair(amp_0, ampc_0);
        amplitude_map[channel + "_paral"] = make_pair(amp_paral, ampc_paral);
        amplitude_map[channel + "_perp"] = make_pair(amp_perp, ampc_perp);
    }
    else if (channel == "Bsdspdsm")
    {
        // b → c(c̄s), spectator s
        amp = lam_bs_c * (getPar("E1t_sccs_BDDb") + getPar("A2t_cscs_BDbD")) - lam_bs_u * (getPar("G1t_scs_BDDb") + getPar("G3t_css_BDDb"));
        amp /= lam_bs_c; // Normalize by lam_bs_c
        ampc = lamst_bs_c * (getPar("E1t_sccs_BDDb") + getPar("A2t_cscs_BDbD")) - lamst_bs_u * (getPar("G1t_scs_BDDb") + getPar("G3t_css_BDDb"));
        ampc /= lamst_bs_c; // Normalize by lamst_bs_c
        amplitude_map[channel] = make_pair(amp, ampc);
    }
    else if (channel == "Bsdpdsm")
    {
        // b → c(c̄d), spectator s
        amp = lam_bd_c * (getPar("E1t_dccs_BDDb")) - lam_bs_u * (getPar("G1t_dcs_BDDb"));
        amp /= lam_bs_c; // Normalize by lam_bs_c
        ampc = lamst_bd_c * (getPar("E1t_dccs_BDDb")) - lamst_bs_u * (getPar("G1t_dcs_BDDb"));
        ampc /= lamst_bs_c; // Normalize by lamst_bs_c
        amplitude_map[channel] = make_pair(amp, ampc);
    }
    else if (channel == "Bsdpdm")
    {
        // b → c(c̄s), spectator s
        amp = lam_bs_c * (getPar("A2t_cdcs_BDbD")) - lam_bs_u * (getPar("G3t_cds_BDDb"));
        amp /= lam_bs_c; // Normalize by lam_bs_c
        ampc = lamst_bs_c * (getPar("A2t_cdcs_BDbD")) - lamst_bs_u * (getPar("G3t_cds_BDDb"));
        ampc /= lamst_bs_c; // Normalize by lamst_bs_c
        amplitude_map[channel] = make_pair(amp, ampc);
    }
    else if (channel == "Bsd0d0b")
    {
        // b → c(c̄s), spectator s
        amp = -lam_bs_c * (getPar("A2t_cdcs_BDbD") + getPar("dP3EW_ucs_BDbD")) -
              lam_bs_u * (getPar("A2_dcds_BDDb") + getPar("dP3EW_ucs_BDbD") - getPar("G3t_cds_BDDb"));
        amp /= lam_bs_c; // Normalize by lam_bs_c
        ampc = -lamst_bs_c * (getPar("A2t_cdcs_BDbD") + getPar("dP3EW_ucs_BDbD")) -
               lamst_bs_u * (getPar("A2_dcds_BDDb") + getPar("dP3EW_ucs_BDbD") - getPar("G3t_cds_BDDb"));
        ampc /= lamst_bs_c; // Normalize by lamst_bs_c
        amplitude_map[channel] = make_pair(amp, ampc);
    }
    else if (channel == "Bddspdsm")
    {
        // b → c(c̄d), spectator d
        amp = lam_bd_c * (getPar("A2t_cscd_BDbD")) - lam_bd_u * (getPar("G3t_csd_BDDb"));
        amp /= lam_bd_c; // Normalize by lam_bd_c
        ampc = lamst_bd_c * (getPar("A2t_cscd_BDbD")) - lamst_bd_u * (getPar("G3t_csd_BDDb"));
        ampc /= lamst_bd_c; // Normalize by lamst_bd_c
        amplitude_map[channel] = make_pair(amp, ampc);
    }
    else if (channel == "Bddspdm")
    {
        // b → c(c̄s), spectator d
        amp = lam_bs_c * (getPar("E1t_sccd_BDDb")) - lam_bs_u * (getPar("G1t_scd_BDDb"));
        amp /= lam_bd_c; // Normalize by lam_bd_c
        ampc = lamst_bs_c * (getPar("E1t_sccd_BDDb")) - lamst_bs_u * (getPar("G1t_scd_BDDb"));
        ampc /= lamst_bd_c; // Normalize by lamst_bd_c
        amplitude_map[channel] = make_pair(amp, ampc);
    }
    else if (channel == "Bddpdm")
    {
        // b → c(c̄d), spectator d
        amp = lam_bd_c * (getPar("E1t_dccd_BDDb") + getPar("A2t_cdcd_BDbD")) -
              lam_bd_u * (getPar("G1t_dcd_BDDb") + getPar("G3t_cdd_BDDb"));
        amp /= lam_bd_c; // Normalize by lam_bd_c
        ampc = lamst_bd_c * (getPar("E1t_dccd_BDDb") + getPar("A2t_cdcd_BDbD")) -
               lamst_bd_u * (getPar("G1t_dcd_BDDb") + getPar("G3t_cdd_BDDb"));
        ampc /= lamst_bd_c; // Normalize by lamst_bd_c
        amplitude_map[channel] = make_pair(amp, ampc);
    }
    else if (channel == "Bdd0d0b")
    {
        // b → c(c̄d), spectator d
        amp = -lam_bd_c * (getPar("A2t_cdcd_BDbD") + getPar("dP3EW_ucd_BDbD")) -
              lam_bd_u * (getPar("A2_dcdd_BDDb") + getPar("dP3EW_ucd_BDbD") - getPar("G3t_cdd_BDDb"));
        ampc = -lamst_bd_c * (getPar("A2t_cdcd_BDbD") + getPar("dP3EW_ucd_BDbD")) -
               lamst_bd_u * (getPar("A2_dcdd_BDDb") + getPar("dP3EW_ucd_BDbD") - getPar("G3t_cdd_BDDb"));
        amp /= lam_bd_c;    // Normalize by lam_bd_c
        ampc /= lamst_bd_c; // Normalize by lamst_bd_c
        amplitude_map[channel] = make_pair(amp, ampc);
    }
    else if (channel == "Bpdpd0b")
    {
        // b → c(c̄d), spectator u
        amp = +lam_bd_c * (getPar("E1t_dccd_BDDb") + getPar("dP1EW_dcu_BDDb")) +
              lam_bd_u * (getPar("A1_dcdd_BDDb") + getPar("dP1EW_dcu_BDDb") - getPar("G1t_dcd_BDDb"));
        ampc = +lamst_bd_c * (getPar("E1t_dccd_BDDb") + getPar("dP1EW_dcu_BDDb")) +
               lamst_bd_u * (getPar("A1_dcdd_BDDb") + getPar("dP1EW_dcu_BDDb") - getPar("G1t_dcd_BDDb"));
        amplitude_map[channel] = make_pair(amp, ampc);
    }
    else if (channel == "Bpdspd0b")
    {
        // b → c(c̄s), spectator u
        amp = +lam_bs_c * (getPar("E1t_sccd_BDDb") + getPar("dP1EW_scu_BDDb")) +
              lam_bs_u * (getPar("A1_scdd_BDDb") + getPar("dP1EW_scu_BDDb") - getPar("G1t_scd_BDDb"));
        ampc = +lamst_bs_c * (getPar("E1t_sccd_BDDb") + getPar("dP1EW_scu_BDDb")) +
               lamst_bs_u * (getPar("A1_scdd_BDDb") + getPar("dP1EW_scu_BDDb") - getPar("G1t_scd_BDDb"));
        amplitude_map[channel] = make_pair(amp, ampc);
    }
    else
    {
        cout << "WARNING: amplitude for channel " << channel << " not found in compute_decay_amplitudes" << endl;
    }
}

// Destructor
goldenmodesB::~goldenmodesB()
{
    // Ensure any dynamically allocated resources are properly cleaned
    // delete histos;  // Uncomment if histos is dynamically allocated
}

//-----------------------------------------------------------
// Helper function to parse a decay channel name
pair<string, pair<string, string>>
goldenmodesB::parseChannel(const string &channel) const
{
    // Identify the B meson part
    string bMeson;
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
        throw runtime_error("Error in parseChannel: Unknown B meson in channel: " + channel);
    }

    string remaining = channel.substr(bMeson.length());

    // Check if it's a J/psi channel or a DD channel
    const string jpsiIdentifier = "jpsi";
    if (remaining.rfind(jpsiIdentifier, 0) == 0)
    {
        // J/psi channel: B→J/ψ X
        // Extract the second final-state meson
        remaining = remaining.substr(jpsiIdentifier.length());
        string meson2;

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
            throw runtime_error("Error in parseChannel: Unable to parse second final-state meson in channel: " + channel);
        }

        return {bMeson, {jpsiIdentifier, meson2}};
    }
    else
    {
        // DD channel: B→D X or B→Ds X
        // Parse the two D mesons from the remaining string
        string meson1, meson2;

        // Try to match known D meson patterns (order matters - check longer patterns first)
        if (remaining.rfind("dspdsm", 0) == 0)
        {
            meson1 = "dsp";
            meson2 = "dsm";
        }
        else if (remaining.rfind("dspd0b", 0) == 0)
        {
            meson1 = "dsp";
            meson2 = "d0b";
        }
        else if (remaining.rfind("dspdm", 0) == 0)
        {
            meson1 = "dsp";
            meson2 = "dm";
        }
        else if (remaining.rfind("dpdsm", 0) == 0)
        {
            meson1 = "dp";
            meson2 = "dsm";
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
        else if (remaining.rfind("d0d0b", 0) == 0)
        {
            meson1 = "d0";
            meson2 = "d0b";
        }
        else
        {
            throw runtime_error("Error in parseChannel: Unable to parse DD channel: " + channel);
        }

        return {bMeson, {meson1, meson2}};
    }
}

//----------------------------------------------------------
// Getter for B meson lifetime
double goldenmodesB::getBMesonLifetime(const string &bMeson) const
{
    static const unordered_map<string, double> lifetimes = {
        {"Bp", tau_Bp},
        {"Bd", tau_Bd},
        {"Bs", tau_Bs}};

    auto it = lifetimes.find(bMeson);
    if (it != lifetimes.end())
    {
        return it->second;
    }

    throw runtime_error("Error in getBMesonLifetime: Unknown B meson '" + bMeson + "'");
}

//----------------------------------------------------------
// Getter for B meson mass
double goldenmodesB::getBMesonMass(const string &bMeson) const
{
    static const unordered_map<string, double> masses = {
        {"Bp", m_Bp},
        {"Bd", m_Bd},
        {"Bs", m_Bs}};

    auto it = masses.find(bMeson);
    if (it != masses.end())
    {
        return it->second;
    }

    throw runtime_error("Error in getBMesonMass: Unknown B meson '" + bMeson + "'");
}

double goldenmodesB::CalculateBR(TComplex amplitude, TComplex amplitude_conj, const string &channel) const
{
    // Parse the channel to extract meson components
    auto parsed = parseChannel(channel);
    const string &bMeson = parsed.first;
    const string &meson1 = parsed.second.first;
    const string &meson2 = parsed.second.second;

    // Get the masses of the decaying B meson and final-state mesons
    double m_B = getBMesonMass(bMeson);
    double m1 = mesonMasses.at(meson1);
    double m2 = mesonMasses.at(meson2);

    // Ensure physical kinematics: prevent sqrt of a negative number
    double mass_term1 = (m_B * m_B - (m1 + m2) * (m1 + m2));
    double mass_term2 = (m_B * m_B - (m1 - m2) * (m1 - m2));

    if (mass_term1 < 0 || mass_term2 < 0)
    {
        throw runtime_error("Error in CalculateBR: Kinematically forbidden decay for channel " + channel);
    }

    // Compute the magnitude of the final-state momentum
    double p = sqrt(mass_term1 * mass_term2) / (2.0 * m_B);

    // add the CKM factor we divided out in the amplitude calculation
    if (bMeson == "Bd")
    {
        amplitude *= lam_bd_c;        // Use lam_bd_c for Bd decays
        amplitude_conj *= lamst_bd_c; // Use lamst_bd_c for conjugate amplitude
    }
    else if (bMeson == "Bs")
    {
        amplitude *= lam_bs_c;        // Use lam_bs_c for Bs decays
        amplitude_conj *= lamst_bs_c; // Use lamst_bs_c for conjugate amplitude
    }

    // Compute the CP-averaged decay width (Γ)
    double decay_width = (G_F * G_F / (32.0 * M_PI * h_t)) * (p / (m_B * m_B)) * (0.5 * (amplitude.Rho2() + amplitude_conj.Rho2()));

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
double goldenmodesB::CalculateAcp(const TComplex &amplitude, const TComplex &conjugate_amplitude) const
{
    // Ensure amplitudes are non-zero to avoid division errors
    double A2 = amplitude.Rho2();
    double Abar2 = conjugate_amplitude.Rho2();

    if (A2 + Abar2 == 0)
    {
        cerr << "Warning: CalculateAcp - Both amplitudes are zero, returning 0." << endl;
        return 0.0;
    }

    // Compute A_CP
    return (Abar2 - A2) / (A2 + Abar2);
}

// ---------------------------------------------------------------------
// Function to calculate direct CP violation parameter C
double goldenmodesB::CalculateC(const TComplex &amplitude, const TComplex &conjugate_amplitude, const string &channel)
{
    // Parse the channel to determine the B meson type
    auto parsed = parseChannel(channel);
    string bMeson = parsed.first;

    // Having factored out the right CKM element in the amplitude calculation, here we just have 2beta or 2betas plus possible NP phase

    TComplex q_p = (bMeson == "Bd") ? TComplex::Exp(TComplex(0, -getParameterValue("myphid"))) : TComplex::Exp(TComplex(0, getParameterValue("myphis")));

    // Special case for K0s and K0l channels: apply q/p_KS, taking into account that the minus sign for the KL is compensated by the CP eigenvalue
    if (channel == "Bdjpsik0s" || channel == "Bdjpsik0l")
    {
        q_p = q_p * ckm.get_q_p_KS(); // Multiply by q/p for K0 mixing
    }
    else if (channel == "Bsjpsik0s" || channel == "Bsjpsik0l")
    {
        q_p = q_p / ckm.get_q_p_KS(); // Multiply by p/q for K0 mixing
    }

    // // Get CP eigenvalue for the channel (ensure it exists)
    // if (cpEigenvalue.find(channel) == cpEigenvalue.end())
    // {
    //     throw runtime_error("Error: CalculateC - CP eigenvalue missing for channel: " + channel);
    // }
    // double eta = cpEigenvalue.at(channel);

    // Compute λ (lambda) = η * (q/p) * (A_conjugate / A)
    if (amplitude.Rho() == 0)
    {
        cerr << "Error: CalculateC - Zero amplitude for " << channel << ", division by zero detected." << endl;
        return 0.0; // Return safe value instead of crashing
    }

    TComplex lambda = q_p * (conjugate_amplitude / amplitude);

    // Compute C observable: C = (1 - |λ|^2) / (1 + |λ|^2)
    double mod_lambda_squared = lambda.Rho2();
    return (1.0 - mod_lambda_squared) / (1.0 + mod_lambda_squared);
}

// ---------------------------------------------------------------------
// Function to calculate CP violation parameter S and Delta S
pair<double, double> goldenmodesB::CalculateS(const TComplex &amplitude, const TComplex &conjugate_amplitude, const string &channel)
{
    // Parse the channel to determine the B meson type
    auto parsed = parseChannel(channel);
    string bMeson = parsed.first;

    bool isBd = (bMeson == "Bd");

    // Having factored out the right CKM element in the amplitude calculation, here we just have 2beta or 2betas plus possible NP phase
    TComplex q_p = isBd ? TComplex::Exp(TComplex(0, -getParameterValue("myphid"))) : TComplex::Exp(TComplex(0, getParameterValue("myphis")));

    // Special case for K0s and K0l channels (apply q/p_KS)
    if (channel == "Bdjpsik0s" || channel == "Bdjpsik0l")
    {
        q_p = q_p * ckm.get_q_p_KS(); // Multiply by q/p for K0 mixing
    }
    else if (channel == "Bsjpsik0s" || channel == "Bsjpsik0l")
    {
        q_p = q_p / ckm.get_q_p_KS(); // Multiply by p/q for K0 mixing
    }

    // // Get CP eigenvalue for the channel (ensure it exists)
    // if (cpEigenvalue.find(channel) == cpEigenvalue.end())
    // {
    //     throw runtime_error("Error: CalculateS - CP eigenvalue missing for channel: " + channel);
    // }
    // double eta = cpEigenvalue.at(channel);

    // Compute λ (lambda) = η * (q/p) * (A_conjugate / A)
    if (amplitude.Rho() == 0)
    {
        cerr << "Error: CalculateS - Zero amplitude for " << channel << ", division by zero detected." << endl;
        return make_pair(0., 0.); // Return safe value instead of crashing
    }

    TComplex lambda = q_p * (conjugate_amplitude / amplitude);

    // Compute S observable: S = 2 Im(λ) / (1 + |λ|^2)
    double mod_lambda_squared = lambda.Rho2();
    double S = -(2.0 * lambda.Im()) / (1.0 + mod_lambda_squared);
    double DeltaS = isBd ? S - sin(getParameterValue("myphid")) : S + sin(getParameterValue("myphis"));

    return make_pair(S, DeltaS);
}

tuple<double, double, double> goldenmodesB::CalculatePhiAndLambda(const TComplex &amplitude, const TComplex &conjugate_amplitude, const string &channel)
{
    // Ensure the amplitude is nonzero to avoid division by zero
    if (amplitude.Rho() == 0)
    {
        throw runtime_error("Error: CalculatePhiAndLambda - Zero amplitude detected for channel: " + channel);
    }

    // Parse the channel to determine the B meson type
    auto parsed = parseChannel(channel);
    string bMeson = parsed.first;

    bool isBd = (bMeson == "Bd");

    // Having factored out the right CKM element in the amplitude calculation, here we just have 2beta or 2betas plus possible NP phase.
    TComplex q_p = isBd ? TComplex::Exp(TComplex(0, -getParameterValue("myphid"))) : TComplex::Exp(TComplex(0, getParameterValue("myphis")));
    double sign = -1.;

    // Compute lambda = (q/p) * (A_cp / A_conj)
    TComplex lambda = q_p * (conjugate_amplitude / amplitude);

    // Compute |lambda|
    double mod_lambda = lambda.Rho();

    // Compute phi = sign * arg(lambda)
    double phi = sign * lambda.Theta();

    // cout << "Channel: " << channel << ", reference angle: " << (isBd ? 2.*ckm.get_beta() : -2.*ckm.get_betas()) << ", calculated phi: " << phi << ", |lambda|: " << mod_lambda << endl;

    return make_tuple(phi, mod_lambda, phi - (isBd ? getParameterValue("myphid") : getParameterValue("myphis")));
}

// pair<vector<string>, string> goldenmodesB::extractChannelFromCorrKey(const string &corr_key)
// {
//     vector<string> channels;
//     string experiment;

//     if (corr_key.rfind("CS_", 0) == 0)
//     {
//         // Format: "CS_channel_exp"
//         size_t underscore = corr_key.find("_", 3);
//         if (underscore != string::npos)
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
//         if (first_underscore != string::npos && second_underscore != string::npos)
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
//         if (underscore != string::npos)
//         {
//             channels.push_back(corr_key.substr(11, underscore - 11));
//             experiment = corr_key.substr(underscore + 1);
//         }
//     }
//     else if (corr_key.rfind("phi_", 0) == 0)
//     {
//         // Format: "phi_channel_exp"
//         size_t underscore = corr_key.find("_", 4);
//         if (underscore != string::npos)
//         {
//             channels.push_back(corr_key.substr(4, underscore - 4));
//             experiment = corr_key.substr(underscore + 1);
//         }
//     }
//     else if (corr_key.rfind("polarization_", 0) == 0)
//     {
//         // Format: "polarization_channel_exp"
//         size_t underscore = corr_key.find("_", 13);
//         if (underscore != string::npos)
//         {
//             channels.push_back(corr_key.substr(13, underscore - 13));
//             experiment = corr_key.substr(underscore + 1);
//         }
//     }
//     else
//     {
//         cerr << "Warning!! Unknown key format: " << corr_key << endl;
//     }

//     return {channels, experiment};
// }

map<string, double> goldenmodesB::CalculatePolarizations(
    pair<TComplex, TComplex> &amplitude_0, pair<TComplex, TComplex> &amplitude_paral, pair<TComplex, TComplex> &amplitude_perp)
{
    map<string, double> polarization_pars;

    // Get decay amplitudes
    TComplex amp_0 = amplitude_0.first;
    TComplex amp_paral = amplitude_paral.first;
    TComplex amp_perp = amplitude_perp.first;
    TComplex conj_amp_0 = amplitude_0.second;
    TComplex conj_amp_paral = amplitude_paral.second;
    TComplex conj_amp_perp = amplitude_perp.second;

    // Calculate squared amplitudes
    double abs2_A0 = amp_0.Rho2();
    double abs2_Aparal = amp_paral.Rho2();
    double abs2_Aperp = amp_perp.Rho2();
    double abs2_conjA0 = conj_amp_0.Rho2();
    double abs2_conjAparal = conj_amp_paral.Rho2();
    double abs2_conjAperp = conj_amp_perp.Rho2();

    // CP-averaged norms
    double norm_A0 = 0.5 * (abs2_A0 + abs2_conjA0);
    double norm_Aparal = 0.5 * (abs2_Aparal + abs2_conjAparal);
    double norm_Aperp = 0.5 * (abs2_Aperp + abs2_conjAperp);

    double norm_amp = norm_A0 + norm_Aparal + norm_Aperp;

    // Calculate CP-averaged polarization fractions
    double f_0 = norm_A0 / norm_amp;
    double f_perp = norm_Aperp / norm_amp;
    double f_paral = norm_Aparal / norm_amp;

    // Calculate CP-averaged relative phases
    double delta_0 = remainder((amp_0.Theta() + conj_amp_0.Theta()) / 2., 2. * M_PI); // Reference phase for A0
    double delta_paral = remainder((amp_paral.Theta() + conj_amp_paral.Theta()) / 2., 2. * M_PI);
    double delta_perp = remainder((amp_perp.Theta() + conj_amp_perp.Theta()) / 2., 2. * M_PI);

    // Store the results in the map
    polarization_pars["f_0"] = f_0;
    polarization_pars["f_perp"] = f_perp;
    polarization_pars["f_paral"] = f_paral;
    polarization_pars["delta_paral"] = delta_paral;
    polarization_pars["delta_perp"] = delta_perp;
    polarization_pars["delta_0"] = delta_0;

    return polarization_pars;
}

double goldenmodesB::Calculate_UncorrelatedObservables(map<string, pair<TComplex, TComplex>> &amplitude_map)
{
    double ll_uncorr = 0.0; // Initialize log-likelihood contribution

    pair<TComplex, TComplex> amp_pair;
    pair<TComplex, TComplex> amp0_pair;
    pair<TComplex, TComplex> ampparal_pair;
    pair<TComplex, TComplex> ampperp_pair;

    vector<pair<string, pair<string, string>>> br_ratios = {
        {"R_Bpjpsipp_Bpjpsikp", {"Bpjpsipp", "Bpjpsikp"}},
        {"R_Bdjpsiom_Bdjpsirho0", {"Bdjpsiom", "Bdjpsirho0"}},
        {"R_Bdjpsikst_Bdjpsik0", {"Bdjpsikst0", "Bdjpsik0"}},
        {"R_Bdjpsieta_Bsjpsieta", {"Bdjpsieta", "Bsjpsieta"}},
        {"R_Bdjpsietap_Bsjpsietap", {"Bdjpsietap", "Bsjpsietap"}},
        {"R_Bdjpsietap_Bdjpsieta", {"Bdjpsietap", "Bdjpsieta"}},
        {"R_Bsjpsieta_Bsjpsiphi", {"Bsjpsieta", "Bsjpsiphi"}},
        {"R_Bsjpsietap_Bsjpsieta", {"Bsjpsietap", "Bsjpsieta"}},
        {"R_Bsjpsietap_Bsjpsiphi", {"Bsjpsietap", "Bsjpsiphi"}}};

    // **Handle Delta A (`deltaA_`)**
    // Format: {ratioKey, {channel1, channel2}} - using channel names, not ACP keys
    vector<pair<string, pair<string, string>>> deltaA_ratios = {
        {"deltaA_Bpjpsipp_Bpjpsikp", {"Bpjpsipp", "Bpjpsikp"}}};

    // Handle differences of phases among different channels from Dalitz analyses
    vector<pair<string, pair<string, string>>> phase_differences = {
        {"Delta_delta_paral_Bdjpsiom-delta_0_Bdjpsirho0", {"delta_paral_Bdjpsiom", "delta_0_Bdjpsirho0"}},
        {"Delta_delta_perp_Bdjpsiom-delta_perp_Bdjpsirho0", {"delta_perp_Bdjpsiom", "delta_perp_Bdjpsirho0"}},
        {"Delta_delta_0_Bdjpsiom-delta_0_Bdjpsirho0", {"delta_0_Bdjpsiom", "delta_0_Bdjpsirho0"}}};

    string basechannel;
    regex pattern(R"(B[dsp][^_-]*)"); // regexp to match the base channel (e.g. Bd0d0b) in the observable name, assuming it starts with B followed by d, s or p and then any characters until an underscore
    smatch match;                  // Object to hold the result
    for (const auto meas_pair : meas)
    {
        const string &key = meas_pair.first;

        // Check if the observable is already computed
        if (obs.find(key) == obs.end())
        {

            if (regex_search(key, match, pattern))
            {
                basechannel = match.str();
                if (Debug)
                {
                    cout << "Extracted base channel: " << basechannel << " from observable key: " << key << endl;
                }
            }
            else
            {
                cerr << "couldn't extract basechannel from " << key << endl;
                continue;
            }

            // Check if it's an averaged vector meson channel
            bool is_vector_channel = find(vectorMesonChannels.begin(), vectorMesonChannels.end(), basechannel) != vectorMesonChannels.end();
            bool is_polarized_measurement = key.find("_0") != string::npos || key.find("_paral") != string::npos || key.find("_perp") != string::npos;

            if (is_vector_channel)
            {
                auto it0 = amplitude_map.find(basechannel + "_0");
                auto itparal = amplitude_map.find(basechannel + "_paral");
                auto itperp = amplitude_map.find(basechannel + "_perp");
                if (it0 == amplitude_map.end() || itparal == amplitude_map.end() || itperp == amplitude_map.end())
                {
                    cerr << "Warning: Polarized amplitudes not found for " << basechannel << endl;
                    continue;
                }
                amp0_pair = it0->second;
                ampparal_pair = itparal->second;
                ampperp_pair = itperp->second;
            }
            else
            {
                auto it = amplitude_map.find(basechannel);
                if (it == amplitude_map.end())
                {
                    cerr << "Warning: Amplitude not found for channel " << basechannel << " in Calculate_CorrelatedObservables" << endl;
                    continue;
                }
                else
                {
                    amp_pair = it->second;
                }
            }

            // Compute the observable based on its type
            if (key.rfind("C", 0) == 0)
            {
                double c_value;
                if (is_vector_channel && !is_polarized_measurement)
                {
                    // Average over polarizations for C
                    double c_0 = CalculateC(amp0_pair.first, amp0_pair.second, basechannel);
                    double c_paral = CalculateC(ampparal_pair.first, ampparal_pair.second, basechannel);
                    double c_perp = CalculateC(ampperp_pair.first, ampperp_pair.second, basechannel);
                    auto polfracs = CalculatePolarizations(amp0_pair, ampparal_pair, ampperp_pair);
                    c_value = c_0 * polfracs["f_0"] + c_paral * polfracs["f_paral"] + c_perp * polfracs["f_perp"];
                }
                else if (is_vector_channel && is_polarized_measurement)
                {
                    // Polarized measurement for vector meson channel
                    if (key.find("_0") != string::npos)
                    {
                        c_value = CalculateC(amp0_pair.first, amp0_pair.second, basechannel);
                    }
                    else if (key.find("_paral") != string::npos)
                    {
                        c_value = CalculateC(ampparal_pair.first, ampparal_pair.second, basechannel);
                    }
                    else if (key.find("_perp") != string::npos)
                    {
                        c_value = CalculateC(ampperp_pair.first, ampperp_pair.second, basechannel);
                    }
                    else
                    {
                        cerr << "Error: Unknown polarization in " << key << endl;
                        continue;
                    }
                }
                else
                {
                    c_value = CalculateC(amp_pair.first, amp_pair.second, basechannel);
                }
                obs[key] = c_value;
            }

            else if (key.rfind("S", 0) == 0)
            {
                double s_value, delta_S;
                if (is_vector_channel && !is_polarized_measurement)
                {
                    // Average over polarizations for S
                    auto s_0 = CalculateS(amp0_pair.first, amp0_pair.second, basechannel);
                    auto s_paral = CalculateS(ampparal_pair.first, ampparal_pair.second, basechannel);
                    auto s_perp = CalculateS(ampperp_pair.first, ampperp_pair.second, basechannel);
                    auto polfracs = CalculatePolarizations(amp0_pair, ampparal_pair, ampperp_pair);
                    s_value = s_0.first * polfracs["f_0"] + s_paral.first * polfracs["f_paral"] + s_perp.first * polfracs["f_perp"];
                    delta_S = s_0.second * polfracs["f_0"] + s_paral.second * polfracs["f_paral"] + s_perp.second * polfracs["f_perp"];
                }
                else if (is_vector_channel && is_polarized_measurement)
                {
                    // Polarized measurement for vector meson channel
                    if (key.find("_0") != string::npos)
                    {
                        auto s_pair = CalculateS(amp0_pair.first, amp0_pair.second, basechannel);
                        s_value = s_pair.first;
                        delta_S = s_pair.second;
                    }
                    else if (key.find("_paral") != string::npos)
                    {
                        auto s_pair = CalculateS(ampparal_pair.first, ampparal_pair.second, basechannel);
                        s_value = s_pair.first;
                        delta_S = s_pair.second;
                    }
                    else if (key.find("_perp") != string::npos)
                    {
                        auto s_pair = CalculateS(ampperp_pair.first, ampperp_pair.second, basechannel);
                        s_value = s_pair.first;
                        delta_S = s_pair.second;
                    }
                    else
                    {
                        cerr << "Error: Unknown polarization in " << key << endl;
                        continue;
                    }
                }
                else
                {
                    auto s_pair = CalculateS(amp_pair.first, amp_pair.second, basechannel);
                    s_value = s_pair.first;
                    delta_S = s_pair.second;
                }
                obs[key] = s_value;
                obs["DeltaS_" + key.substr(1)] = delta_S;
            }
            else if (key.rfind("phis", 0) == 0)
            {
                if (is_vector_channel && !is_polarized_measurement)
                {
                    // Average over polarizations for phi_s
                    auto phi_lambda_0 = CalculatePhiAndLambda(amp0_pair.first, amp0_pair.second, basechannel);
                    auto phi_lambda_paral = CalculatePhiAndLambda(ampparal_pair.first, ampparal_pair.second, basechannel);
                    auto phi_lambda_perp = CalculatePhiAndLambda(ampperp_pair.first, ampperp_pair.second, basechannel);
                    auto polfracs = CalculatePolarizations(amp0_pair, ampparal_pair, ampperp_pair);

                    double phi_0 = get<0>(phi_lambda_0);
                    double phi_paral = get<0>(phi_lambda_paral);
                    double phi_perp = get<0>(phi_lambda_perp);

                    double phi_avg = phi_0 * polfracs["f_0"] + phi_paral * polfracs["f_paral"] + phi_perp * polfracs["f_perp"];
                    double delta_phi_avg = get<2>(phi_lambda_0) * polfracs["f_0"] + get<2>(phi_lambda_paral) * polfracs["f_paral"] + get<2>(phi_lambda_perp) * polfracs["f_perp"];
                    obs[key] = phi_avg;
                    obs["Delta" + key] = delta_phi_avg;
                }
                else if (is_vector_channel && is_polarized_measurement)
                {
                    // Polarized measurement for vector meson channel
                    double phi_pol, delta_phi_pol;
                    if (key.find("_0") != string::npos)
                    {
                        auto phi_lambda_0 = CalculatePhiAndLambda(amp0_pair.first, amp0_pair.second, basechannel);
                        phi_pol = get<0>(phi_lambda_0);
                        delta_phi_pol = get<2>(phi_lambda_0);
                    }
                    else if (key.find("_paral") != string::npos)
                    {
                        auto phi_lambda_paral = CalculatePhiAndLambda(ampparal_pair.first, ampparal_pair.second, basechannel);
                        phi_pol = get<0>(phi_lambda_paral);
                        delta_phi_pol = get<2>(phi_lambda_paral);
                    }
                    else if (key.find("_perp") != string::npos)
                    {
                        auto phi_lambda_perp = CalculatePhiAndLambda(ampperp_pair.first, ampperp_pair.second, basechannel);
                        phi_pol = get<0>(phi_lambda_perp);
                        delta_phi_pol = get<2>(phi_lambda_perp);
                    }
                    else
                    {
                        cerr << "Error: Unknown polarization in " << key << endl;
                        continue;
                    }
                    obs[key] = phi_pol;
                    obs["Delta" + key] = delta_phi_pol;
                }
                else
                {
                    // Non-vector meson channel
                    auto phi_lambda_result = CalculatePhiAndLambda(amp_pair.first, amp_pair.second, basechannel);
                    double phi = get<0>(phi_lambda_result);
                    double delta_phi = get<2>(phi_lambda_result);
                    obs[key] = phi;
                    obs["Delta" + key] = delta_phi;
                }
            }
            else if (key.rfind("lambda", 0) == 0)
            {
                if (is_vector_channel && !is_polarized_measurement)
                {
                    // Average over polarizations for |lambda|
                    auto phi_lambda_0 = CalculatePhiAndLambda(amp0_pair.first, amp0_pair.second, basechannel);
                    auto phi_lambda_paral = CalculatePhiAndLambda(ampparal_pair.first, ampparal_pair.second, basechannel);
                    auto phi_lambda_perp = CalculatePhiAndLambda(ampperp_pair.first, ampperp_pair.second, basechannel);
                    auto polfracs = CalculatePolarizations(amp0_pair, ampparal_pair, ampperp_pair);

                    double lambda_0 = get<1>(phi_lambda_0);
                    double lambda_paral = get<1>(phi_lambda_paral);
                    double lambda_perp = get<1>(phi_lambda_perp);

                    double lambda_avg = lambda_0 * polfracs["f_0"] + lambda_paral * polfracs["f_paral"] + lambda_perp * polfracs["f_perp"];
                    obs[key] = lambda_avg;
                }
                else if (is_vector_channel && is_polarized_measurement)
                {
                    // Polarized measurement for vector meson channel
                    double lambda_pol;
                    if (key.find("_0") != string::npos)
                    {
                        auto phi_lambda_0 = CalculatePhiAndLambda(amp0_pair.first, amp0_pair.second, basechannel);
                        lambda_pol = get<1>(phi_lambda_0);
                    }
                    else if (key.find("_paral") != string::npos)
                    {
                        auto phi_lambda_paral = CalculatePhiAndLambda(ampparal_pair.first, ampparal_pair.second, basechannel);
                        lambda_pol = get<1>(phi_lambda_paral);
                    }
                    else if (key.find("_perp") != string::npos)
                    {
                        auto phi_lambda_perp = CalculatePhiAndLambda(ampperp_pair.first, ampperp_pair.second, basechannel);
                        lambda_pol = get<1>(phi_lambda_perp);
                    }
                    else
                    {
                        cerr << "Error: Unknown polarization in " << key << endl;
                        continue;
                    }
                    obs[key] = lambda_pol;
                }
                else
                {
                    // Non-vector meson channel
                    auto phi_lambda_result = CalculatePhiAndLambda(amp_pair.first, amp_pair.second, basechannel);
                    double lambda_val = get<1>(phi_lambda_result);
                    obs[key] = lambda_val;
                }
            }
            else if (key.rfind("f_", 0) == 0 || key.rfind("delta_", 0) == 0)
            {
                // Compute polarization parameters for vector meson channels
                auto pol_params = CalculatePolarizations(amp0_pair, ampparal_pair, ampperp_pair);
                obs["f_0_" + basechannel] = pol_params.at("f_0");
                obs["f_paral_" + basechannel] = pol_params.at("f_paral");
                obs["f_perp_" + basechannel] = pol_params.at("f_perp");
                obs["delta_paral_" + basechannel] = pol_params.at("delta_paral")-pol_params.at("delta_0");
                obs["delta_perp_" + basechannel] = pol_params.at("delta_perp")-pol_params.at("delta_0");
            }
            else if (key.rfind("2beta", 0) == 0)
            {
                if (is_vector_channel && !is_polarized_measurement)
                {
                    cerr << "Error: 2beta observable not implemented for unpolarized vector meson channels: " << basechannel << endl;
                }
                else if (is_vector_channel && is_polarized_measurement)
                {
                    // Polarized measurement for vector meson channel
                    double twobeta_pol, delta_twobeta_pol;
                    if (key.find("_0") != string::npos)
                    {
                        auto phi_lambda_0 = CalculatePhiAndLambda(amp0_pair.first, amp0_pair.second, basechannel);
                        twobeta_pol = get<0>(phi_lambda_0);
                        delta_twobeta_pol = get<2>(phi_lambda_0);
                    }
                    else if (key.find("_paral") != string::npos)
                    {
                        auto phi_lambda_paral = CalculatePhiAndLambda(ampparal_pair.first, ampparal_pair.second, basechannel);
                        auto phi_lambda_0 = CalculatePhiAndLambda(amp0_pair.first, amp0_pair.second, basechannel);
                        twobeta_pol = remainder(get<0>(phi_lambda_paral) - get<0>(phi_lambda_0), 2. * M_PI);
                        delta_twobeta_pol = remainder(get<2>(phi_lambda_paral) + get<2>(phi_lambda_0), 2. * M_PI);
                    }
                    else if (key.find("_perp") != string::npos)
                    {
                        auto phi_lambda_perp = CalculatePhiAndLambda(ampperp_pair.first, ampperp_pair.second, basechannel);
                        auto phi_lambda_0 = CalculatePhiAndLambda(amp0_pair.first, amp0_pair.second, basechannel);
                        twobeta_pol = remainder(get<0>(phi_lambda_perp) - get<0>(phi_lambda_0), 2. * M_PI);
                        delta_twobeta_pol = remainder(get<2>(phi_lambda_perp) + get<2>(phi_lambda_0), 2. * M_PI);
                    }
                    else
                    {
                        cerr << "Error: Unknown polarization in " << key << endl;
                        continue;
                    }
                    obs[key] = twobeta_pol;
                    obs["Delta" + key] = delta_twobeta_pol;
                }
                else
                {
                    // Non-vector meson channel
                    auto phi_lambda_result = CalculatePhiAndLambda(amp_pair.first, amp_pair.second, basechannel);
                    double phi_s = get<0>(phi_lambda_result);
                    double delta_phi_s = get<2>(phi_lambda_result);

                    obs[key] = phi_s;
                    obs["Delta" + key] = delta_phi_s;
                }
            }
            else if (key.rfind("alpha_", 0) == 0)
            {
                if (is_vector_channel && !is_polarized_measurement)
                {
                    cerr << "Error: alpha observable not implemented for unpolarized vector meson channels: " << basechannel << endl;
                }
                else if (is_vector_channel && is_polarized_measurement)
                {
                    // Polarized measurement for vector meson channel
                    double alpha_pol;
                    if (key.find("0") != string::npos)
                    {
                        double lambda = get<1>(CalculatePhiAndLambda(amp0_pair.first, amp0_pair.second, basechannel));
                        alpha_pol = (1. - lambda) / (1. + lambda);
                    }
                    else if (key.find("paral") != string::npos)
                    {
                        double lambda = get<1>(CalculatePhiAndLambda(ampparal_pair.first, ampparal_pair.second, basechannel));
                        alpha_pol = (1. - lambda) / (1. + lambda);
                    }
                    else if (key.find("perp") != string::npos)
                    {
                        double lambda = get<1>(CalculatePhiAndLambda(ampperp_pair.first, ampperp_pair.second, basechannel));
                        alpha_pol = (1. - lambda) / (1. + lambda);
                    }
                    else
                    {
                        cerr << "Error: Unknown polarization in " << key << endl;
                        continue;
                    }
                    obs[key] = alpha_pol;
                }
                else
                {
                    cerr << "Error: alpha observable not implemented for non-vector meson channels: " << basechannel << endl;
                }
            }
            else if (key.rfind("ACP", 0) == 0)
            {
                double acp_value;
                if (is_vector_channel && !is_polarized_measurement)
                {
                    // Average over polarizations for ACP
                    double acp_0 = CalculateAcp(amp0_pair.first, amp0_pair.second);
                    double acp_paral = CalculateAcp(ampparal_pair.first, ampparal_pair.second);
                    double acp_perp = CalculateAcp(ampperp_pair.first, ampperp_pair.second);
                    auto polfracs = CalculatePolarizations(amp0_pair, ampparal_pair, ampperp_pair);
                    acp_value = acp_0 * polfracs["f_0"] + acp_paral * polfracs["f_paral"] + acp_perp * polfracs["f_perp"];
                }
                else if (is_vector_channel && is_polarized_measurement)
                {
                    // Polarized measurement for vector meson channel
                    if (key.find("0") != string::npos)
                    {
                        acp_value = CalculateAcp(amp0_pair.first, amp0_pair.second);
                    }
                    else if (key.find("paral") != string::npos)
                    {
                        acp_value = CalculateAcp(ampparal_pair.first, ampparal_pair.second);
                    }
                    else if (key.find("perp") != string::npos)
                    {
                        acp_value = CalculateAcp(ampperp_pair.first, ampperp_pair.second);
                    }
                    else
                    {
                        cerr << "Error: Unknown polarization in " << key << endl;
                        continue;
                    }
                }
                else
                {
                    acp_value = CalculateAcp(amp_pair.first, amp_pair.second);
                }
                obs[key] = acp_value;
            }
            else if (key.find("BR") == 0)
            {
                double br_value = 0.;
                if (is_vector_channel)
                {
                    // Average over polarizations for BR
                    br_value += CalculateBR(amp0_pair.first, amp0_pair.second, basechannel);
                    br_value += CalculateBR(ampparal_pair.first, ampparal_pair.second, basechannel);
                    br_value += CalculateBR(ampperp_pair.first, ampperp_pair.second, basechannel);
                }
                else
                {
                    br_value = CalculateBR(amp_pair.first, amp_pair.second, basechannel);
                }
                obs[key] = br_value;
            }
            else if (key.find("R_") == 0)
            {
                // Check if this is one of the defined BR ratios
                auto ratio_it = find_if(br_ratios.begin(), br_ratios.end(), [&key](const pair<string, pair<string, string>> &ratio)
                                        { return key == ratio.first; });

                if (ratio_it != br_ratios.end())
                {
                    const string &numerator_channel = ratio_it->second.first;
                    const string &denominator_channel = ratio_it->second.second;

                    double br_numerator, br_denominator;

                    is_vector_channel = find(vectorMesonChannels.begin(), vectorMesonChannels.end(), numerator_channel) != vectorMesonChannels.end();
                    is_polarized_measurement = numerator_channel.find("_0") != string::npos || numerator_channel.find("_paral") != string::npos || numerator_channel.find("_perp") != string::npos;

                    if (is_vector_channel)
                    {
                        amp0_pair = amplitude_map.at(numerator_channel + "_0");
                        ampparal_pair = amplitude_map.at(numerator_channel + "_paral");
                        ampperp_pair = amplitude_map.at(numerator_channel + "_perp");
                    }
                    else
                    {
                        amp_pair = amplitude_map.at(numerator_channel);
                    }
                    // Calculate BR for numerator
                    if (is_vector_channel && !is_polarized_measurement)
                    {
                        br_numerator = CalculateBR(amp0_pair.first, amp0_pair.second, numerator_channel) +
                                       CalculateBR(ampparal_pair.first, ampparal_pair.second, numerator_channel) +
                                       CalculateBR(ampperp_pair.first, ampperp_pair.second, numerator_channel);
                    }
                    else if (is_vector_channel && is_polarized_measurement)
                    {
                        if (numerator_channel.find("_0") != string::npos)
                        {
                            br_numerator = CalculateBR(amp0_pair.first, amp0_pair.second, numerator_channel);
                        }
                        else if (numerator_channel.find("_paral") != string::npos)
                        {
                            br_numerator = CalculateBR(ampparal_pair.first, ampparal_pair.second, numerator_channel);
                        }
                        else if (numerator_channel.find("_perp") != string::npos)
                        {
                            br_numerator = CalculateBR(ampperp_pair.first, ampperp_pair.second, numerator_channel);
                        }
                        else
                        {
                            cerr << "Error: Unknown polarization in numerator channel " << numerator_channel << " for BR ratio " << key << endl;
                            continue;
                        }
                    }
                    else
                    {
                        br_numerator = CalculateBR(amp_pair.first, amp_pair.second, numerator_channel);
                    }

                    // Calculate BR for denominator

                    is_vector_channel = find(vectorMesonChannels.begin(), vectorMesonChannels.end(), denominator_channel) != vectorMesonChannels.end();
                    is_polarized_measurement = denominator_channel.find("_0") != string::npos || denominator_channel.find("_paral") != string::npos || denominator_channel.find("_perp") != string::npos;

                    if (is_vector_channel)
                    {
                        amp0_pair = amplitude_map.at(denominator_channel + "_0");
                        ampparal_pair = amplitude_map.at(denominator_channel + "_paral");
                        ampperp_pair = amplitude_map.at(denominator_channel + "_perp");
                    }
                    else
                    {
                        amp_pair = amplitude_map.at(denominator_channel);
                    }
                    // Calculate BR for denominator
                    if (is_vector_channel && !is_polarized_measurement)
                    {
                        br_denominator = CalculateBR(amp0_pair.first, amp0_pair.second, denominator_channel) +
                                         CalculateBR(ampparal_pair.first, ampparal_pair.second, denominator_channel) +
                                         CalculateBR(ampperp_pair.first, ampperp_pair.second, denominator_channel);
                    }
                    else if (is_vector_channel && is_polarized_measurement)
                    {
                        if (denominator_channel.find("_0") != string::npos)
                        {
                            br_denominator = CalculateBR(amp0_pair.first, amp0_pair.second, denominator_channel);
                        }
                        else if (denominator_channel.find("_paral") != string::npos)
                        {
                            br_denominator = CalculateBR(ampparal_pair.first, ampparal_pair.second, denominator_channel);
                        }
                        else if (denominator_channel.find("_perp") != string::npos)
                        {
                            br_denominator = CalculateBR(ampperp_pair.first, ampperp_pair.second, denominator_channel);
                        }
                        else
                        {
                            cerr << "Error: Unknown polarization in denominator channel " << denominator_channel << " for BR ratio " << key << endl;
                            continue;
                        }
                    }
                    else
                    {
                        br_denominator = CalculateBR(amp_pair.first, amp_pair.second, denominator_channel);
                    }

                    if (br_denominator != 0)
                    {
                        obs[key] = br_numerator / br_denominator;
                    }
                    else
                    {
                        cerr << "Error: Denominator BR is zero for ratio " << key << endl;
                        continue;
                    }
                }
                else
                {
                    cerr << "Error: Unknown BR ratio " << key << endl;
                    continue;
                }
            }
            else if (key.find("deltaA_") == 0)
            {
                // Check if this is one of the defined Delta A ratios
                auto ratio_it = find_if(deltaA_ratios.begin(), deltaA_ratios.end(), [&key](const pair<string, pair<string, string>> &ratio)
                                        { return key == ratio.first; });

                if (ratio_it != deltaA_ratios.end())
                {
                    const string &channel1 = ratio_it->second.first;
                    const string &channel2 = ratio_it->second.second;

                    double acp1, acp2;

                    is_vector_channel = find(vectorMesonChannels.begin(), vectorMesonChannels.end(), channel1) != vectorMesonChannels.end();
                    is_polarized_measurement = channel1.find("_0") != string::npos || channel1.find("_paral") != string::npos || channel1.find("_perp") != string::npos;

                    if (is_vector_channel)
                    {
                        amp0_pair = amplitude_map.at(channel1 + "_0");
                        ampparal_pair = amplitude_map.at(channel1 + "_paral");
                        ampperp_pair = amplitude_map.at(channel1 + "_perp");
                    }
                    else
                    {
                        amp_pair = amplitude_map.at(channel1);
                    }


                    // Calculate ACP for channel 1
                    if (is_vector_channel && !is_polarized_measurement)
                    {
                        double acp_0 = CalculateAcp(amp0_pair.first, amp0_pair.second);
                        double acp_paral = CalculateAcp(ampparal_pair.first, ampparal_pair.second);
                        double acp_perp = CalculateAcp(ampperp_pair.first, ampperp_pair.second);
                        auto polfracs = CalculatePolarizations(amp0_pair, ampparal_pair, ampperp_pair);
                        acp1 = acp_0 * polfracs["f_0"] + acp_paral * polfracs["f_paral"] + acp_perp * polfracs["f_perp"];
                    }
                    else
                    {
                        acp1 = CalculateAcp(amp_pair.first, amp_pair.second);
                    }

                    is_vector_channel = find(vectorMesonChannels.begin(), vectorMesonChannels.end(), channel2) != vectorMesonChannels.end();
                    is_polarized_measurement = channel2.find("_0") != string::npos || channel2.find("_paral") != string::npos || channel2.find("_perp") != string::npos;

                    if (is_vector_channel)
                    {
                        amp0_pair = amplitude_map.at(channel2 + "_0");
                        ampparal_pair = amplitude_map.at(channel2 + "_paral");
                        ampperp_pair = amplitude_map.at(channel2 + "_perp");
                    }
                    else
                    {
                        amp_pair = amplitude_map.at(channel2);
                    }


                    // Calculate ACP for channel 2
                    if (is_vector_channel && !is_polarized_measurement)
                    {
                        double acp_0 = CalculateAcp(amp0_pair.first, amp0_pair.second);
                        double acp_paral = CalculateAcp(ampparal_pair.first, ampparal_pair.second);
                        double acp_perp = CalculateAcp(ampperp_pair.first, ampperp_pair.second);
                        auto polfracs = CalculatePolarizations(amp0_pair, ampparal_pair, ampperp_pair);
                        acp2 = acp_0 * polfracs["f_0"] + acp_paral * polfracs["f_paral"] + acp_perp * polfracs["f_perp"];
                    }
                    else
                    {
                        acp2 = CalculateAcp(amp_pair.first, amp_pair.second);
                    }
                    obs[key] = acp1 - acp2;
                }
                else
                {
                    cerr << "Error: Unknown Delta A ratio " << key << endl;
                    continue;
                }
            }
            else if (key.find("Delta_delta_") == 0)
            {
                // Check if this is one of the defined phase differences
                auto ratio_it = find_if(phase_differences.begin(), phase_differences.end(), [&key](const pair<string, pair<string, string>> &ratio)
                                        { return key == ratio.first; });

                if (ratio_it != phase_differences.end())
                {
                    const string &channel1 = ratio_it->second.first;
                    const string &channel2 = ratio_it->second.second;

                    if (regex_search(channel1, match, pattern))
                    {
                        basechannel = match.str();
                        if (Debug)
                        {
                            cout << "Extracted base channel: " << basechannel << " from observable key: " << key << endl;
                        }
                    }
                    else
                    {
                        cerr << "couldn't extract basechannel from " << key << endl;
                        continue;
                    }

                    string poltype;
                    if (channel1.find("_0") != string::npos)
                    {
                        poltype = "_0";
                    }
                    else if (channel1.find("_paral") != string::npos)
                    {
                        poltype = "_paral";
                    }
                    else if (channel1.find("_perp") != string::npos)
                    {
                        poltype = "_perp";
                    }
                    else
                    {
                        cerr << "Error: Unknown polarization type in " << channel1 << endl;
                        continue;
                    }

                    double delta1, delta2;

                    amp0_pair = amplitude_map.at(basechannel + "_0");
                    ampparal_pair = amplitude_map.at(basechannel + "_paral");
                    ampperp_pair = amplitude_map.at(basechannel + "_perp");


                    auto pols1 = CalculatePolarizations(amp0_pair, ampparal_pair, ampperp_pair);
                    delta1 = pols1.at("delta_" + poltype.substr(1)); // Extract the relevant delta based on polarization type

                    if (regex_search(channel2, match, pattern))
                    {
                        basechannel = match.str();
                        if (Debug)
                        {
                            cout << "Extracted base channel: " << basechannel << " from observable key: " << key << endl;
                        }
                    }
                    else
                    {
                        cerr << "couldn't extract basechannel from " << key << endl;
                        continue;
                    }

                    if (channel2.find("_0") != string::npos)
                    {
                        poltype = "_0";
                    }
                    else if (channel2.find("_paral") != string::npos)
                    {
                        poltype = "_paral";
                    }
                    else if (channel2.find("_perp") != string::npos)
                    {
                        poltype = "_perp";
                    }
                    else
                    {
                        cerr << "Error: Unknown polarization type in " << channel2 << endl;
                        continue;
                    }

                    amp0_pair = amplitude_map.at(basechannel + "_0");
                    ampparal_pair = amplitude_map.at(basechannel + "_paral");
                    ampperp_pair = amplitude_map.at(basechannel + "_perp");
                    auto pols2 = CalculatePolarizations(amp0_pair, ampparal_pair, ampperp_pair);
                    delta2 = pols2.at("delta_" + poltype.substr(1)); // Extract the relevant delta based on polarization type

                    // Calculate the phase difference and ensure it's within [-pi, pi]
                    double delta_diff = remainder(delta1 - delta2, 2. * M_PI);
                    obs[key] = delta_diff;
                }
                else
                {                     cerr << "Error: Unknown phase difference " << key << endl;
                    continue;
                }
            }
            else
            {
                cerr << "Error: Unknown observable " << key << endl;
                continue;
            }
        }
        else
        {
            if (Debug)
            {
                cout << "Observable " << key << " already computed: " << obs[key] << endl;
            }
        }

        ll_uncorr += meas_pair.second.logweight(obs[key]);
        if (Debug)
        {
            cout << "Observable: " << key << " Value: " << obs[key] << " Measured: " << meas_pair.second.getMean() << " Uncertainty: " << meas_pair.second.getSigma() << " isAngle: " << meas_pair.second.getIsAngle() << endl;
        }
    }

    return ll_uncorr;
}

//----------------------------------------------------------------------------------

double goldenmodesB::Calculate_CorrelatedObservables(map<string, pair<TComplex, TComplex>> &amplitude_map)
{
    double ll_corr = 0.0; // Initialize the log-likelihood contribution

    string basechannel;
    for (const auto meas_pair : corrmeas)
    {
        const string &key = meas_pair.first;
        const CorrelatedGaussianObservables &corrObs = meas_pair.second;
        const vector<string> &obs_names = corrmeas_channels.at(key);

        if (obs_names.size() != corrObs.getNObs())
        {
            cerr << "Error: Mismatch in number of observables for " << key << endl;
            continue;
        }

        for (const auto &obs_name : obs_names)
        {
            // Check if the observable is already computed
            if (obs.find(obs_name) == obs.end())
            {

                // Parse the channel from the observable name
                size_t pos = obs_name.find("B");
                if (pos != string::npos)
                {
                    basechannel = obs_name.substr(pos);
                }
                else
                {
                    cerr << "Delimiter not found in " << obs_name << " so couldn't extract basechannel" << endl;
                    continue;
                }

                // Check if it's an averaged vector meson channel
                bool is_vector_channel = find(vectorMesonChannels.begin(), vectorMesonChannels.end(), basechannel) != vectorMesonChannels.end();
                bool is_polarized_measurement = obs_name.find("0") != string::npos || obs_name.find("paral") != string::npos || obs_name.find("perp") != string::npos;

                pair<TComplex, TComplex> amp_pair;
                pair<TComplex, TComplex> amp0_pair;
                pair<TComplex, TComplex> ampparal_pair;
                pair<TComplex, TComplex> ampperp_pair;

                if (is_vector_channel)
                {
                    auto it0 = amplitude_map.find(basechannel + "_0");
                    auto itparal = amplitude_map.find(basechannel + "_paral");
                    auto itperp = amplitude_map.find(basechannel + "_perp");
                    if (it0 == amplitude_map.end() || itparal == amplitude_map.end() || itperp == amplitude_map.end())
                    {
                        cerr << "Warning: Polarized amplitudes not found for " << basechannel << endl;
                        continue;
                    }
                    amp0_pair = it0->second;
                    ampparal_pair = itparal->second;
                    ampperp_pair = itperp->second;
                }
                else
                {
                    auto it = amplitude_map.find(basechannel);
                    if (it == amplitude_map.end())
                    {
                        cerr << "Warning: Amplitude not found for channel " << basechannel << " in Calculate_CorrelatedObservables" << endl;
                        continue;
                    }
                    else
                    {
                        amp_pair = it->second;
                    }
                }

                // Compute the observable based on its type
                if (obs_name.rfind("C", 0) == 0)
                {
                    double c_value;
                    if (is_vector_channel && !is_polarized_measurement)
                    {
                        // Average over polarizations for C
                        double c_0 = CalculateC(amp0_pair.first, amp0_pair.second, basechannel);
                        double c_paral = CalculateC(ampparal_pair.first, ampparal_pair.second, basechannel);
                        double c_perp = CalculateC(ampperp_pair.first, ampperp_pair.second, basechannel);
                        auto polfracs = CalculatePolarizations(amp0_pair, ampparal_pair, ampperp_pair);
                        c_value = c_0 * polfracs["f_0"] + c_paral * polfracs["f_paral"] + c_perp * polfracs["f_perp"];
                    }
                    else if (is_vector_channel && is_polarized_measurement)
                    {
                        // Polarized measurement for vector meson channel
                        if (obs_name.find("0") != string::npos)
                        {
                            c_value = CalculateC(amp0_pair.first, amp0_pair.second, basechannel);
                        }
                        else if (obs_name.find("paral") != string::npos)
                        {
                            c_value = CalculateC(ampparal_pair.first, ampparal_pair.second, basechannel);
                        }
                        else if (obs_name.find("perp") != string::npos)
                        {
                            c_value = CalculateC(ampperp_pair.first, ampperp_pair.second, basechannel);
                        }
                        else
                        {
                            cerr << "Error: Unknown polarization in " << obs_name << endl;
                            continue;
                        }
                    }
                    else
                    {
                        c_value = CalculateC(amp_pair.first, amp_pair.second, basechannel);
                    }
                    obs[obs_name] = c_value;
                }

                else if (obs_name.rfind("S", 0) == 0)
                {
                    double s_value, delta_S;
                    if (is_vector_channel && !is_polarized_measurement)
                    {
                        // Average over polarizations for S
                        auto s_0 = CalculateS(amp0_pair.first, amp0_pair.second, basechannel);
                        auto s_paral = CalculateS(ampparal_pair.first, ampparal_pair.second, basechannel);
                        auto s_perp = CalculateS(ampperp_pair.first, ampperp_pair.second, basechannel);
                        auto polfracs = CalculatePolarizations(amp0_pair, ampparal_pair, ampperp_pair);
                        s_value = s_0.first * polfracs["f_0"] + s_paral.first * polfracs["f_paral"] + s_perp.first * polfracs["f_perp"];
                        delta_S = s_0.second * polfracs["f_0"] + s_paral.second * polfracs["f_paral"] + s_perp.second * polfracs["f_perp"];
                    }
                    else if (is_vector_channel && is_polarized_measurement)
                    {
                        // Polarized measurement for vector meson channel
                        if (obs_name.find("0") != string::npos)
                        {
                            auto s_pair = CalculateS(amp0_pair.first, amp0_pair.second, basechannel);
                            s_value = s_pair.first;
                            delta_S = s_pair.second;
                        }
                        else if (obs_name.find("paral") != string::npos)
                        {
                            auto s_pair = CalculateS(ampparal_pair.first, ampparal_pair.second, basechannel);
                            s_value = s_pair.first;
                            delta_S = s_pair.second;
                        }
                        else if (obs_name.find("perp") != string::npos)
                        {
                            auto s_pair = CalculateS(ampperp_pair.first, ampperp_pair.second, basechannel);
                            s_value = s_pair.first;
                            delta_S = s_pair.second;
                        }
                        else
                        {
                            cerr << "Error: Unknown polarization in " << obs_name << endl;
                            continue;
                        }
                    }
                    else
                    {
                        auto s_pair = CalculateS(amp_pair.first, amp_pair.second, basechannel);
                        s_value = s_pair.first;
                        delta_S = s_pair.second;
                    }
                    obs[obs_name] = s_value;
                    obs["DeltaS_" + obs_name.substr(1)] = delta_S;
                }
                else if (obs_name.rfind("phis", 0) == 0)
                {
                    if (is_vector_channel && !is_polarized_measurement)
                    {
                        // Average over polarizations for phi_s
                        auto phi_lambda_0 = CalculatePhiAndLambda(amp0_pair.first, amp0_pair.second, basechannel);
                        auto phi_lambda_paral = CalculatePhiAndLambda(ampparal_pair.first, ampparal_pair.second, basechannel);
                        auto phi_lambda_perp = CalculatePhiAndLambda(ampperp_pair.first, ampperp_pair.second, basechannel);
                        auto polfracs = CalculatePolarizations(amp0_pair, ampparal_pair, ampperp_pair);

                        double phi_0 = get<0>(phi_lambda_0);
                        double phi_paral = get<0>(phi_lambda_paral);
                        double phi_perp = get<0>(phi_lambda_perp);

                        double phi_avg = phi_0 * polfracs["f_0"] + phi_paral * polfracs["f_paral"] + phi_perp * polfracs["f_perp"];
                        double delta_phi_avg = get<2>(phi_lambda_0) * polfracs["f_0"] + get<2>(phi_lambda_paral) * polfracs["f_paral"] + get<2>(phi_lambda_perp) * polfracs["f_perp"];
                        obs[obs_name] = phi_avg;
                        obs["Delta" + obs_name] = delta_phi_avg;
                    }
                    else if (is_vector_channel && is_polarized_measurement)
                    {
                        // Polarized measurement for vector meson channel
                        double phi_pol, delta_phi_pol;
                        if (obs_name.find("0") != string::npos)
                        {
                            auto phi_lambda_0 = CalculatePhiAndLambda(amp0_pair.first, amp0_pair.second, basechannel);
                            phi_pol = get<0>(phi_lambda_0);
                            delta_phi_pol = get<2>(phi_lambda_0);
                        }
                        else if (obs_name.find("paral") != string::npos)
                        {
                            auto phi_lambda_paral = CalculatePhiAndLambda(ampparal_pair.first, ampparal_pair.second, basechannel);
                            phi_pol = get<0>(phi_lambda_paral);
                            delta_phi_pol = get<2>(phi_lambda_paral);
                        }
                        else if (obs_name.find("perp") != string::npos)
                        {
                            auto phi_lambda_perp = CalculatePhiAndLambda(ampperp_pair.first, ampperp_pair.second, basechannel);
                            phi_pol = get<0>(phi_lambda_perp);
                            delta_phi_pol = get<2>(phi_lambda_perp);
                        }
                        else
                        {
                            cerr << "Error: Unknown polarization in " << obs_name << endl;
                            continue;
                        }
                        obs[obs_name] = phi_pol;
                        obs["Delta" + obs_name] = delta_phi_pol;
                    }
                    else
                    {
                        // Non-vector meson channel
                        auto phi_lambda_result = CalculatePhiAndLambda(amp_pair.first, amp_pair.second, basechannel);
                        double phi = get<0>(phi_lambda_result);
                        double delta_phi = get<2>(phi_lambda_result);
                        obs[obs_name] = phi;
                        obs["Delta" + obs_name] = delta_phi;
                    }
                }
                else if (obs_name.rfind("lambda", 0) == 0)
                {
                    if (is_vector_channel && !is_polarized_measurement)
                    {
                        // Average over polarizations for |lambda|
                        auto phi_lambda_0 = CalculatePhiAndLambda(amp0_pair.first, amp0_pair.second, basechannel);
                        auto phi_lambda_paral = CalculatePhiAndLambda(ampparal_pair.first, ampparal_pair.second, basechannel);
                        auto phi_lambda_perp = CalculatePhiAndLambda(ampperp_pair.first, ampperp_pair.second, basechannel);
                        auto polfracs = CalculatePolarizations(amp0_pair, ampparal_pair, ampperp_pair);

                        double lambda_0 = get<1>(phi_lambda_0);
                        double lambda_paral = get<1>(phi_lambda_paral);
                        double lambda_perp = get<1>(phi_lambda_perp);

                        double lambda_avg = lambda_0 * polfracs["f_0"] + lambda_paral * polfracs["f_paral"] + lambda_perp * polfracs["f_perp"];
                        obs[obs_name] = lambda_avg;
                    }
                    else if (is_vector_channel && is_polarized_measurement)
                    {
                        // Polarized measurement for vector meson channel
                        double lambda_pol;
                        if (obs_name.find("0") != string::npos)
                        {
                            auto phi_lambda_0 = CalculatePhiAndLambda(amp0_pair.first, amp0_pair.second, basechannel);
                            lambda_pol = get<1>(phi_lambda_0);
                        }
                        else if (obs_name.find("paral") != string::npos)
                        {
                            auto phi_lambda_paral = CalculatePhiAndLambda(ampparal_pair.first, ampparal_pair.second, basechannel);
                            lambda_pol = get<1>(phi_lambda_paral);
                        }
                        else if (obs_name.find("perp") != string::npos)
                        {
                            auto phi_lambda_perp = CalculatePhiAndLambda(ampperp_pair.first, ampperp_pair.second, basechannel);
                            lambda_pol = get<1>(phi_lambda_perp);
                        }
                        else
                        {
                            cerr << "Error: Unknown polarization in " << obs_name << endl;
                            continue;
                        }
                        obs[obs_name] = lambda_pol;
                    }
                    else
                    {
                        // Non-vector meson channel
                        auto phi_lambda_result = CalculatePhiAndLambda(amp_pair.first, amp_pair.second, basechannel);
                        double lambda_val = get<1>(phi_lambda_result);
                        obs[obs_name] = lambda_val;
                    }
                }
                else if (obs_name.rfind("f_", 0) == 0 || obs_name.rfind("delta_", 0) == 0)
                {
                    // Compute polarization parameters for vector meson channels
                    auto pol_params = CalculatePolarizations(amp0_pair, ampparal_pair, ampperp_pair);
                    obs["f_0_" + basechannel] = pol_params.at("f_0");
                    obs["f_paral_" + basechannel] = pol_params.at("f_paral");
                    obs["f_perp_" + basechannel] = pol_params.at("f_perp");
                    obs["delta_paral_" + basechannel] = pol_params.at("delta_paral") - pol_params.at("delta_0") ; // Store relative phases
                    obs["delta_perp_" + basechannel] = pol_params.at("delta_perp") - pol_params.at("delta_0") ;
                }
                else if (obs_name.rfind("2beta", 0) == 0)
                {
                    if (is_vector_channel && !is_polarized_measurement)
                    {
                        cerr << "Error: 2beta observable not implemented for unpolarized vector meson channels: " << basechannel << endl;
                    }
                    else if (is_vector_channel && is_polarized_measurement)
                    {
                        // Polarized measurement for vector meson channel
                        double twobeta_pol, delta_twobeta_pol;
                        if (obs_name.find("_0") != string::npos)
                        {
                            auto phi_lambda_0 = CalculatePhiAndLambda(amp0_pair.first, amp0_pair.second, basechannel);
                            twobeta_pol = get<0>(phi_lambda_0);
                            delta_twobeta_pol = get<2>(phi_lambda_0);
                        }
                        else if (obs_name.find("_paral") != string::npos)
                        {
                            auto phi_lambda_paral = CalculatePhiAndLambda(ampparal_pair.first, ampparal_pair.second, basechannel);
                            auto phi_lambda_0 = CalculatePhiAndLambda(amp0_pair.first, amp0_pair.second, basechannel);
                            twobeta_pol = remainder(get<0>(phi_lambda_paral) - get<0>(phi_lambda_0), 2. * M_PI);
                            delta_twobeta_pol = remainder(get<2>(phi_lambda_paral) + get<2>(phi_lambda_0), 2. * M_PI);
                        }
                        else if (obs_name.find("_perp") != string::npos)
                        {
                            auto phi_lambda_perp = CalculatePhiAndLambda(ampperp_pair.first, ampperp_pair.second, basechannel);
                            auto phi_lambda_0 = CalculatePhiAndLambda(amp0_pair.first, amp0_pair.second, basechannel);
                            twobeta_pol = remainder(get<0>(phi_lambda_perp) - get<0>(phi_lambda_0), 2. * M_PI);
                            delta_twobeta_pol = remainder(get<2>(phi_lambda_perp) + get<2>(phi_lambda_0), 2. * M_PI);
                        }
                        else
                        {
                            cerr << "Error: Unknown polarization in " << obs_name << endl;
                            continue;
                        }
                        obs[obs_name] = twobeta_pol;
                        obs["Delta" + obs_name] = delta_twobeta_pol;
                    }
                    else
                    {
                        // Non-vector meson channel
                        auto phi_lambda_result = CalculatePhiAndLambda(amp_pair.first, amp_pair.second, basechannel);
                        double phi_s = get<0>(phi_lambda_result);
                        double delta_phi_s = get<2>(phi_lambda_result);

                        obs[obs_name] = phi_s;
                        obs["Delta" + obs_name] = delta_phi_s;
                    }
                }
                else if (obs_name.rfind("alpha_", 0) == 0)
                {
                    if (is_vector_channel && !is_polarized_measurement)
                    {
                        cerr << "Error: alpha observable not implemented for unpolarized vector meson channels: " << basechannel << endl;
                    }
                    else if (is_vector_channel && is_polarized_measurement)
                    {
                        // Polarized measurement for vector meson channel
                        double alpha_pol;
                        if (obs_name.find("0") != string::npos)
                        {
                            double lambda = get<1>(CalculatePhiAndLambda(amp0_pair.first, amp0_pair.second, basechannel));
                            alpha_pol = (1. - lambda) / (1. + lambda);
                        }
                        else if (obs_name.find("paral") != string::npos)
                        {
                            double lambda = get<1>(CalculatePhiAndLambda(ampparal_pair.first, ampparal_pair.second, basechannel));
                            alpha_pol = (1. - lambda) / (1. + lambda);
                        }
                        else if (obs_name.find("perp") != string::npos)
                        {
                            double lambda = get<1>(CalculatePhiAndLambda(ampperp_pair.first, ampperp_pair.second, basechannel));
                            alpha_pol = (1. - lambda) / (1. + lambda);
                        }
                        else
                        {
                            cerr << "Error: Unknown polarization in " << obs_name << endl;
                            continue;
                        }
                        obs[obs_name] = alpha_pol;
                    }
                    else
                    {
                        cerr << "Error: alpha observable not implemented for non-vector meson channels: " << basechannel << endl;
                    }
                }
                else if (obs_name.rfind("ACP", 0) == 0)
                {
                    double acp_value;
                    if (is_vector_channel && !is_polarized_measurement)
                    {
                        // Average over polarizations for ACP
                        double acp_0 = CalculateAcp(amp0_pair.first, amp0_pair.second);
                        double acp_paral = CalculateAcp(ampparal_pair.first, ampparal_pair.second);
                        double acp_perp = CalculateAcp(ampperp_pair.first, ampperp_pair.second);
                        auto polfracs = CalculatePolarizations(amp0_pair, ampparal_pair, ampperp_pair);
                        acp_value = acp_0 * polfracs["f_0"] + acp_paral * polfracs["f_paral"] + acp_perp * polfracs["f_perp"];
                    }
                    else if (is_vector_channel && is_polarized_measurement)
                    {
                        // Polarized measurement for vector meson channel
                        if (obs_name.find("0") != string::npos)
                        {
                            acp_value = CalculateAcp(amp0_pair.first, amp0_pair.second);
                        }
                        else if (obs_name.find("paral") != string::npos)
                        {
                            acp_value = CalculateAcp(ampparal_pair.first, ampparal_pair.second);
                        }
                        else if (obs_name.find("perp") != string::npos)
                        {
                            acp_value = CalculateAcp(ampperp_pair.first, ampperp_pair.second);
                        }
                        else
                        {
                            cerr << "Error: Unknown polarization in " << obs_name << endl;
                            continue;
                        }
                    }
                    else
                    {
                        acp_value = CalculateAcp(amp_pair.first, amp_pair.second);
                    }
                    obs[obs_name] = acp_value;
                }
                else
                {
                    cerr << "Error: Unknown observable " << obs_name << endl;
                    continue;
                }
            }
        }

        TVectorD res(corrObs.getNObs());
        for (size_t i = 0; i < corrObs.getNObs(); ++i)
        {
            const string &obs_name = obs_names[i];
            res(i) = obs[obs_name];
        }

        if (Debug)
            cout << "computing logweight for correlated gaussian observable " << meas_pair.first << endl;
        ll_corr += corrObs.logweight(res);
    }

    return ll_corr;
}

//------------------------------------------------------------

double goldenmodesB::LogLikelihood(const vector<double> &parameters)
{

    obs.clear();     // Clear obs map for each iteration
    double ll = 0.0; // Log-likelihood accumulator

    if (parameters.empty())
    {
        cerr << "Error: Empty parameters vector!" << endl;
        return log(0.);
    }

    // Unpack CKM parameters and compute CKM elements
    vector<double> ckmParams(parameters.end() - 5, parameters.end() - 1);
    ckm.computeCKM(ckmParams[0], ckmParams[1], ckmParams[2], ckmParams[3], true);

    obs["2beta"] = ckm.get_beta() * 2.0;
    obs["phis"] = ckm.get_betas() * -2.0;

    // eta-eta' mixing angle
    double theta_P = parameters.back();
    SetParameterValue("theta_P", theta_P);

    lam_bs_c = ckm.getVcs() * TComplex::Conjugate(ckm.getVcb());
    lam_bs_u = ckm.getVus() * TComplex::Conjugate(ckm.getVub());
    lam_bd_c = ckm.getVcd() * TComplex::Conjugate(ckm.getVcb());
    lam_bd_u = ckm.getVud() * TComplex::Conjugate(ckm.getVub());
    lamst_bs_c = TComplex::Conjugate(lam_bs_c);
    lamst_bs_u = TComplex::Conjugate(lam_bs_u);
    lamst_bd_c = TComplex::Conjugate(lam_bd_c);
    lamst_bd_u = TComplex::Conjugate(lam_bd_u);

    // Populate parameterValues from BAT parameters directly
    // This ensures proper mapping regardless of channelParameters ordering
    for (unsigned int i = 0; i < GetNParameters(); ++i)
    {
        const string &paramName = GetParameter(i).GetName();
        SetParameterValue(paramName, parameters[i]);
        // assign to obs for histogram filling
        obs[paramName] = parameters[i];
    }
    for (const auto &channel : channelNamesSU3)
    {
        // fill histos for mod and phase for the effective parameters
        for (const auto &param : channelParameters[channel])
        {

            string newStr;
            size_t length = param.length();

            newStr.reserve(length + 1);

            if (length >= 3 && param.substr(length - 3) == "_re")
            {
                newStr.append(param, 0, length - 3); // Append text before "_re"
                // strip the possible trailing "delta_" from the parameter name to get the full parameter
                newStr.erase(0, newStr.find("delta_") == 0 ? 6 : 0); // Remove "delta_" if it exists at the start
                obs[newStr + "_abs"] = getPar(newStr).Rho();         // Append replacement
                obs[newStr + "_arg"] = getPar(newStr).Theta();       // Append replacement
            }
        }
    }

    // Build the amplitude map using physical amplitudes

    amplitude_map.clear();

    for (const string &channel : channelNamesSU3)
    {
        compute_decay_amplitudes(channel);

        // For CP observables, add K0s and K0l cases with the same amplitudes as Bdjpsik0
        if (channel == "Bdjpsik0")
        {
            amplitude_map["Bdjpsik0s"] = make_pair(amplitude_map.at("Bdjpsik0").first / sqrt(2.), amplitude_map.at("Bdjpsik0").second / sqrt(2.));
            amplitude_map["Bdjpsik0l"] = make_pair(amplitude_map.at("Bdjpsik0").first / sqrt(2.), amplitude_map.at("Bdjpsik0").second / sqrt(2.));
            //    cout << "Mapped Bdjpsik0 to Bdjpsik0s and Bdjpsik0l" << endl;
        }
        else if (channel == "Bsjpsik0b")
        {
            // Also map Bsjpsik0b -> Bsjpsik0s for observables
            amplitude_map["Bsjpsik0s"] = make_pair(amplitude_map.at("Bsjpsik0b").first / sqrt(2.), amplitude_map.at("Bsjpsik0b").second / sqrt(2.));
            amplitude_map["Bsjpsik0l"] = make_pair(-amplitude_map.at("Bsjpsik0b").first / sqrt(2.), -amplitude_map.at("Bsjpsik0b").second / sqrt(2.));
            //    cout << "Mapped Bsjpsik0b to Bsjpsik0s" << endl;
        }
    }

    for (const string &channel : channels)
    {
        // compute amplitudes for eta and etaprime final states
        if ((channel.length() >= strlen("eta") && channel.compare(channel.length() - strlen("eta"), strlen("eta"), "eta") == 0))
        {
            amplitude_map[channel] = make_pair(cos(theta_P) * amplitude_map.at(channel + "8").first - sin(theta_P) * amplitude_map.at(channel + "1").first,
                                               cos(theta_P) * amplitude_map.at(channel + "8").second - sin(theta_P) * amplitude_map.at(channel + "1").second);
            amplitude_map[channel + "p"] = make_pair(sin(theta_P) * amplitude_map.at(channel + "8").first + cos(theta_P) * amplitude_map.at(channel + "1").first,
                                                     sin(theta_P) * amplitude_map.at(channel + "8").second + cos(theta_P) * amplitude_map.at(channel + "1").second);
        }
    }

    // Check NaN/Inf values in amplitudes
    for (const auto &[chan, amp_pair] : amplitude_map)
    {
        if (isnan(amp_pair.first.Rho()) || isinf(amp_pair.first.Rho()) ||
            isnan(amp_pair.second.Rho()) || isinf(amp_pair.second.Rho()))
        {
            cerr << "Invalid amplitude (NaN or Inf) for channel: " << chan << endl;
            return log(0.);
        }
    }

    // Calculate BRs, polarization fractions, CP asymmetries, etc., and fill obs map
    for (const string &channel : channels)
    {
        auto ch = parseChannel(channel);

        bool is_vector_channel = find(vectorMesonChannels.begin(), vectorMesonChannels.end(), channel) != vectorMesonChannels.end();

        // Branching Ratio
        double br = CalculateBR(amplitude_map[channel].first, amplitude_map[channel].second, channel);
        obs["BR_" + channel] = br;
        // Polarization fractions and phases for vector meson channels
        if (is_vector_channel)
        {
            auto pol_params = CalculatePolarizations(amplitude_map[channel + "_0"], amplitude_map[channel + "_paral"], amplitude_map[channel + "_perp"]);
            obs["f_0_" + channel] = pol_params.at("f_0");
            obs["f_paral_" + channel] = pol_params.at("f_paral");
            obs["f_perp_" + channel] = pol_params.at("f_perp");
            obs["delta_paral_" + channel] = pol_params.at("delta_paral") - pol_params.at("delta_0"); // Store relative phases
            obs["delta_perp_" + channel] = pol_params.at("delta_perp") - pol_params.at("delta_0");
        }
        // Direct CP Asymmetry for charged B decays
        if (ch.first == "Bp")
        {
            if (is_vector_channel)
            {
                // Average over polarizations for ACP in vector meson channels
                double acp_0 = CalculateAcp(amplitude_map[channel + "_0"].first, amplitude_map[channel + "_0"].second);
                obs["ACP_" + channel + "_0"] = acp_0;
                double acp_paral = CalculateAcp(amplitude_map[channel + "_paral"].first, amplitude_map[channel + "_paral"].second);
                obs["ACP_" + channel + "_paral"] = acp_paral;
                double acp_perp = CalculateAcp(amplitude_map[channel + "_perp"].first, amplitude_map[channel + "_perp"].second);
                obs["ACP_" + channel + "_perp"] = acp_perp;
                double acp_avg = acp_0 * obs["f_0_" + channel] + acp_paral * obs["f_paral_" + channel] + acp_perp * obs["f_perp_" + channel];
                obs["ACP_" + channel] = acp_avg;
            }
            else
            {
                double acp = CalculateAcp(amplitude_map[channel].first, amplitude_map[channel].second);
                obs["ACP_" + channel] = acp;
            }
        }
        // S and C parameters for neutral B decays
        if (ch.first == "Bd" || ch.first == "Bs")
        {
            if (is_vector_channel)
            {
                // Average over polarizations for S and C in vector meson channels
                double c_0 = CalculateC(amplitude_map[channel + "_0"].first, amplitude_map[channel + "_0"].second, channel);
                auto s_0 = CalculateS(amplitude_map[channel + "_0"].first, amplitude_map[channel + "_0"].second, channel);
                double c_paral = CalculateC(amplitude_map[channel + "_paral"].first, amplitude_map[channel + "_paral"].second, channel);
                auto s_paral = CalculateS(amplitude_map[channel + "_paral"].first, amplitude_map[channel + "_paral"].second, channel);
                double c_perp = CalculateC(amplitude_map[channel + "_perp"].first, amplitude_map[channel + "_perp"].second, channel);
                auto s_perp = CalculateS(amplitude_map[channel + "_perp"].first, amplitude_map[channel + "_perp"].second, channel);

                double c_avg = c_0 * obs["f_0_" + channel] + c_paral * obs["f_paral_" + channel] + c_perp * obs["f_perp_" + channel];
                double s_avg = s_0.first * obs["f_0_" + channel] + s_paral.first * obs["f_paral_" + channel] + s_perp.first * obs["f_perp_" + channel];
                double delta_S_avg = s_0.second * obs["f_0_" + channel] + s_paral.second * obs["f_paral_" + channel] + s_perp.second * obs["f_perp_" + channel];

                obs["C_" + channel] = c_avg;
                obs["S_" + channel] = s_avg;
                obs["DeltaS_" + channel] = delta_S_avg;
            }
            else
            {
                double c_param = CalculateC(amplitude_map[channel].first, amplitude_map[channel].second, channel);
                auto s_pair = CalculateS(amplitude_map[channel].first, amplitude_map[channel].second, channel);
                obs["C_" + channel] = c_param;
                obs["S_" + channel] = s_pair.first;
                obs["DeltaS_" + channel] = s_pair.second;
            }

            // auto phi_lambda = CalculatePhiAndLambda(amplitude_map[channel].first, amplitude_map[channel].second, channel);
            // obs["phi_" + channel] = get<0>(phi_lambda);
            // obs["lambda_" + channel] = get<1>(phi_lambda);
            // obs["Deltaphi_" + channel] = get<2>(phi_lambda);
        }
    }

    // Add contributions from uncorrelated and correlated observables
    ll += Calculate_UncorrelatedObservables(amplitude_map);
    ll += Calculate_CorrelatedObservables(amplitude_map);

    if (Debug)
    {
        cout << "LogLikelihood: " << ll << endl;
        for (const auto &[obs_name, obs_value] : obs)
        {
            cout << obs_name << ": " << obs_value << endl;
        }
    }

    return ll;
}

//---------------------------------------------------------

void goldenmodesB::MCMCUserIterationInterface()
{

    vector<double> pars;

    for (unsigned int i = 0; i < fMCMCNChains; ++i)
    {
        pars = fMCMCStates.at(i).parameters;
        try
        {
            obs["LogLikelihood"] = LogLikelihood(pars);
        }
        catch (const exception &e)
        {
            cerr << "Error in LogLikelihood: " << e.what() << endl;
            return;
        }

        // Directly fill histograms
        histos.fillh1d();
        histos.fillh2d();
    }
}

void goldenmodesB::SaveHistograms(const string &filename)
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
    cout << "Histograms saved to " << filename << endl;
}

// ---------------------------------------------------------

// void goldenmodesB::PrintObservablePulls(const string &filename)
// {
//     ofstream outfile(filename);

//     if (!outfile.is_open())
//     {
//         cerr << "Error: Unable to open file " << filename << endl;
//         return;
//     }

//     // Header for the output file
//     outfile << "Observable\t\t Measurement\t\t Mean\t\t Sigma\t\t Pull\n";

//     // Define mappings for observable names and measurement keys
//     map<string, string> measMap = {
//         {"BR_Bdjpsik0", "BRBdjpsik0"}, {"BR_Bdjpsip0", "BRBdjpsip0"}, {"BR_Bdjpsiom", "BRBdjpsiom"}, {"BR_Bpjpsikp", "BRBpjpsikp"}, {"BR_Bpjpsipp", "BRBpjpsipp"}, {"BR_Bsjpsiph", "BRBsjpsiph"}, {"BR_Bsjpsik0s", "BRBsjpsik0s"}, {"BR_Bdjpsirh", "BRBdjpsirh"}, {"BR_Bdjpsikst", "BRBdjpsikst"}, {"BR_Bsjpsikst", "BRBsjpsikst"}, {"C_Bdjpsik0s", "CBdjpsik0s"}, {"C_Bdjpsik0l", "CBdjpsik0l"}, {"S_Bdjpsik0s", "SBdjpsik0s"}, {"S_Bdjpsik0l", "SBdjpsik0l"}, {"C_Bdjpsip0", "CBdjpsip0"}, {"S_Bdjpsip0", "SBdjpsip0"}, {"ACP_Bpjpsikp", "ACPBpjpsikp"}, {"ACP_Bpjpsipp", "ACPBpjpsipp"}, {"deltaA_Bpjpsipp_Bpjpsikp", "deltaA_Bpjpsipp_Bpjpsikp"}, {"R_Bpjpsipp_Bpjpsikp", "R_Bpjpsipp_Bpjpsikp"}, {"lambda_Bsjpsiph", "lambda_Bsjpsiph"}, {"phi_Bsjpsiph", "phi_Bsjpsiph"}, {"f_perp_Bsjpsiph", "f_perp_Bsjpsiph"}, {"f_paral_Bsjpsiph", "f_paral_Bsjpsiph"}, {"f_0_Bsjpsiph", "f_0_Bsjpsiph"}, {"delta_perp_Bsjpsiph", "delta_perp_Bsjpsiph"}, {"delta_paral_Bsjpsiph", "delta_paral_Bsjpsiph"}, {"R_Bdjpsiom_Bdjpsirh", "R_Bdjpsiom_Bdjpsirh"}, {"R_Bdjpsikst_Bdjpsik0", "R_Bdjpsikst_Bdjpsik0"}, {"C_Bdjpsikst", "CBdjpsikst"}, {"S_Bdjpsikst", "SBdjpsikst"}, {"f_0_Bdjpsikst", "f_0_Bdjpsikst"}, {"f_paral_Bdjpsikst", "f_paral_Bdjpsikst"}, {"f_perp_Bdjpsikst", "f_perp_Bdjpsikst"}, {"delta_paral_Bdjpsikst", "delta_paral_Bdjpsikst"}, {"delta_perp_Bdjpsikst", "delta_perp_Bdjpsikst"}, {"f_0_Bsjpsikst", "f_0_Bsjpsikst"}, {"f_paral_Bsjpsikst", "f_paral_Bsjpsikst"}, {"delta_paral_Bsjpsikst", "delta_paral_Bsjpsikst"}, {"delta_perp_Bsjpsikst", "delta_perp_Bsjpsikst"}, {"A0_CP_Bsjpsikst", "A0_CP_Bsjpsikst"}, {"Aparal_CP_Bsjpsikst", "Aparal_CP_Bsjpsikst"}, {"Aperp_CP_Bsjpsikst", "Aperp_CP_Bsjpsikst"}, {"BR_Bddpdm", "BRBddpdm"}, {"C_Bddpdm", "CBddpdm"}, {"S_Bddpdm", "SBddpdm"}, {"BR_Bsdspdsm", "BRBsdspdsm"}, {"C_Bsdspdsm", "CBsdspdsm"}, {"S_Bsdspdsm", "SBsdspdsm"}, {"BR_Bpdpd0b", "BRBpdpd0b"}, {"ACP_Bpdpd0b", "ACPBpdpd0b"}, {"BR_Bpdspd0b", "BRBpdspd0b"}, {"ACP_Bpdspd0b", "ACPBpdspd0b"}

//     };

//     for (const auto &histPair : histos.h1d)
//     {
//         const string &obs_name = histPair.first;
//         TH1D *hist = histPair.second;

//         double obs_mean = hist->GetMean();
//         double obs_measurement = 0.0;
//         double sigma_measurement = 0.0;
//         bool found_measurement = false;

//         // Check if the observable exists in `newmeas` first
//         if (measMap.find(obs_name) != measMap.end())
//         {
//             const string &key = measMap[obs_name];

//             if (newmeas.find(key) != newmeas.end())
//             {
//                 obs_measurement = newmeas.at(key).getMean();
//                 sigma_measurement = newmeas.at(key).getSigma();
//                 found_measurement = true;
//             }
//             else if (meas.find(key) != meas.end())
//             {
//                 // Use `meas` only if it is not found in `newmeas`
//                 obs_measurement = meas.at(key).getMean();
//                 sigma_measurement = meas.at(key).getSigma();
//                 found_measurement = true;
//             }
//         }

//         if (found_measurement)
//         {
//             // Calculate the pull value: (predicted - observed) / uncertainty
//             double pull = (obs_mean - obs_measurement) / sigma_measurement;

//             outfile << obs_name << "\t"
//                     << obs_measurement << "\t"
//                     << obs_mean << "\t"
//                     << sigma_measurement << "\t"
//                     << pull << "\n";
//         }
//     }

//     outfile.close();
//     cout << "Pull values saved to " << filename << endl;
// }
