#include "goldenmodes.h"
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

goldenmodes::goldenmodes() : BCModel() , histos(obs)
{
  cout << "constructor for goldenmodes called: inserting the experimental data" << endl;
    std::vector<dato> data;
    PDGAverage pdgaverage;


    DeclareParameters();  // Ensure parameters are defined


    // Add CKM parameters directly in the constructor
    AddParameter("CKM_Vud", 0.97415, 0.97447);         // Vud parameter
    AddParameter("CKM_Vcb", 0.04046, 0.04194 );    // Vcb parameter
    AddParameter("CKM_Vub", 0.00349 , 0.00419);  // Vub parameter
    AddParameter("CKM_gamma", 1.12224671, 1.22347581);   // gamma parameter (in rad)


    //Add mixing angle between om1 and om8
    //AddParameter("theta_om", 0.6370451769779, 0.6370451769779); //in rad

    SetPriorConstantAll();

    //measurments

    //BRBdjpsik0
    data.push_back(dato(9.02e-4, 0.10e-4, 0.26e-4)); //Belle:2019xld
    data.push_back(dato(8.1e-4, 0.9e-4, 0.6e-4)); //Belle:2019avj
    data.push_back(dato(8.85e-4, 1.35e-4, 0.1e-4)); //BaBar:2007esv
    data.push_back(dato(8.69e-4, 0.22e-4, 0.30e-4)); //BaBar:2004htr
    data.push_back(dato(9.5e-4, 0.8e-4, 0.6e-4)); //CLEO:2000emb
    data.push_back(dato(11.5e-4, 2.3e-4, 1.7e-4)); //CDF:1995izg
    data.push_back(dato(6.93e-4, 4.07e-4, 0.04e-4)); //CLEO:1991roe
    data.push_back(dato(9.24e-4, 7.21e-4, 0.05e-4)); //ARGUS:1990jet


    pdgaverage.setData(data);
    pdgaverage.setName("BRBdjpsik0");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();


    //BRBdjpsip0
    data.push_back(dato(1.62e-5, 0.11e-5, 0.06e-5)); //Belle:2018nxw
    data.push_back(dato(1.69e-5, 0.14e-5, 0.07e-5)); //BaBar:2008kfx
    data.push_back(dato(2.6e-5, 1.0e-5, 0.2e-5)); //CLEO:2000emb


    pdgaverage.setData(data);
    pdgaverage.setName("BRBdjpsip0");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();


    //BRBdjpsiom
    data.push_back(dato(1.9e-5, 0.6e-5, 0.1e-5)); //LHCb:2014vbo
    data.push_back(dato(2.16e-5, 0.30e-5, 0.14e-5)); //Belle-II:2024fgp

    pdgaverage.setData(data);
    pdgaverage.setName("BRBdjpsiom");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    //BRBpjpsikp
    data.push_back(dato(10.32e-3, 0.07e-3, 0.24e-3));  //BELLE:2019xld
    data.push_back(dato(9.4e-3, 0.7e-3, 0.8e-3)); //Belle:2019avj
    data.push_back(dato(8.9e-3, 0.6e-3, 0.5e-3)); //Belle:2017psv
    data.push_back(dato(8.1e-3, 1.3e-3, 0.7e-3)); //BaBar:2005pcw
    data.push_back(dato(10.61e-3, 0.15e-3, 0.48e-3)); //BaBar:2004htr
    data.push_back(dato(10.4e-3, 1.1e-3, 0.1e-3)); //BaBar:2005sdl
    data.push_back(dato(10.2e-3, 0.8e-3, 0.7e-3)); //CLEO:1997ilq
    data.push_back(dato(9.24e-3, 3.04e-3, 0.05e-3)); //CLEO:1991roe
    data.push_back(dato(8.09e-3, 3.50e-3, 0.04e-3)); //ARGUS:1990jet

    pdgaverage.setData(data);
    pdgaverage.setName("BRBpjpsikp");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    //BRBpjpsipp
    data.push_back(dato(3.8e-5, 0.6e-5, 0.3e-5));  //Belle:2002oex
    data.push_back(dato(2.02e-5, 0.12e-5, 0.10e-5)); //Belle-II:2024hqw


    pdgaverage.setData(data);
    pdgaverage.setName("BRBpjpsipp");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();


    //BRBsjpsiph
    data.push_back(dato(1.037e-3, 0.032e-3, 0.022e-3));  //LHCb:2021qbv
    data.push_back(dato(1.25e-3, 0.07e-3, 0.23e-3)); //Belle:2013sdi
    data.push_back(dato(1.5e-3, 0.5e-3, 0.1e-3)); //CDF:1996ivk


    pdgaverage.setData(data);
    pdgaverage.setName("BRBsjpsiph");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    //BRBdjpsikst
    data.push_back(dato(1.19e-3, 0.01e-3, 0.08e-3));  //Belle:2014nuw
    data.push_back(dato(1.335e-3, 0.215e-3, 0.02e-3)); //BaBar:2007esv
    data.push_back(dato(1.309e-3, 0.026e-3, 0.077e-3)); //BaBar:2004htr
    data.push_back(dato(1.29e-3, 0.05e-3, 0.13e-3)); //Belle:2002otd
    data.push_back(dato(1.74e-3, 0.20e-3, 0.18e-3)); //CDF:1998tqc
    data.push_back(dato(1.32e-3, 0.17e-3, 0.17e-3)); //CLEO:1997ilq
    data.push_back(dato(1.27e-3, 0.65e-3, 0.01e-3)); //CLEO:1991roe
    data.push_back(dato(1.27e-3, 0.60e-3, 0.01e-3)); //ARGUS:1990jet
    data.push_back(dato(4.04e-3, 1.81e-3, 0.02e-3)); //CLEO:1987iba

    pdgaverage.setData(data);
    pdgaverage.setName("BRBdjpsikst");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    //BRBdjpsirh
    data.push_back(dato(2.515e-5, 0.10e-5, 0.165e-5)); //LHCb:2014vbo
    data.push_back(dato(2.7e-5, 0.3e-5, 0.2e-5)); //BaBar:2007yvx

    pdgaverage.setData(data);
    pdgaverage.setName("BRBdjpsirh");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    //BRBsjpsikst
    meas.insert(pair<string, dato>("BRBsjpsikst", dato(4.14e-5, 0.18e-5, 0.35e-5)));

    //BR Ratios

    //R_Bpjpsipp_Bpjpsikp
    data.push_back(dato(3.846e-2, 0.018e-2, 0.018e-2));  //LHCb:2024exp
    data.push_back(dato(3.5e-2, 0.3e-2, 1.2e-2)); //ATLAS:2016rxw
    data.push_back(dato(4.86e-2, 0.82e-2, 0.15e-2)); //CDF:2007mkw
    data.push_back(dato(5.37e-2, 0.45e-2, 0.11e-2)); //BaBar:2004kla
    data.push_back(dato(5.1e-2, 1.8e-2, 0.1e-2));  //CDF:1996efe
    data.push_back(dato(5.2e-2, 2.4e-2)); //CLEO:1995zgs


    pdgaverage.setData(data);
    pdgaverage.setName("R_Bpjpsipp_Bpjpsikp");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();


    data.push_back(dato(0.86, 0.19, 0.10));  //LHCb:2012cw
    data.push_back(dato(0.70, 0.30, 0.05)); //LHCb:2013dkk


    pdgaverage.setData(data);
    pdgaverage.setName("R_Bdjpsiom_Bdjpsirh");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    //R_Bdjpsikst_Bdjpsik0

    data.push_back(dato(1.51, 0.05, 0.08));  //BaBar:2004htr
    data.push_back(dato(1.39, 0.36, 0.10)); //CDF:1996ivk
    pdgaverage.setData(data);

    pdgaverage.setName("R_Bdjpsikst_Bdjpsik0");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    //R_BRBsjpsik0s_Bdjpsik0s
    meas.insert(pair<string, dato>( "R_BRBsjpsik0s_Bdjpsik0s", dato(0.0431, 0.0017, 0.0012, 0.0025));  //LHCb:2015brj

    //ACP measurments

    //CBdjpsip0

    data.push_back(dato(0.155, 0.14, 0.035));  ////Belle:2018nxw
    data.push_back(dato(0.13, 0.12, 0.03)); //Belle-II:2024hqw


    pdgaverage.setData(data);
    pdgaverage.setName("CBdjpsip0");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    //SBdjpsip0
    data.push_back(dato(-0.59, 0.19, 0.03));  ////Belle:2018nxw
    data.push_back(dato(-0.88, 0.17, 0.03)); //Belle-II:2024hqw


    pdgaverage.setData(data);
    pdgaverage.setName("SBdjpsip0");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    meas.insert(pair<string, dato>("CBdjpsik0s", dato(0.005, 0.021, 0.037))); //Belle:2012paq


    //SBdjpsik0s
    data.push_back(dato(0.670, 0.029, 0.013));  ////Belle:2012paq
    data.push_back(dato(0.57, 0.58, 0.06)); //Belle:2012teq
    data.push_back(dato(0.775, 0.425)); //CDF:1999ijp
    data.push_back(dato(0.73, 0.93, 0.16)); //ALEPH:2000jem
    data.push_back(dato(3.1, 1.9, 0.5)); //OPAL:1998mtm


    pdgaverage.setData(data);
    pdgaverage.setName("SBdjpsik0s");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();



    meas.insert(pair<string, dato>("CBdjpsik0l", dato(-0.031, 0.026, 0.029))); //Belle:2012paq
    meas.insert(pair<string, dato>("SBdjpsik0l", dato(0.642, 0.047, 0.021))); //Belle:2012paq


    //SBsjpsik0s
    meas.insert(pair<string, dato>("SBsjpsik0s", dato(-0.08, 0.40, 0.08))); //LHCb:2015brj
    meas.insert(pair<string, dato>("CBsjpsik0s", dato(-0.28, 0.41, 0.08))); //LHCb:2015brj


    //ACPBpjpsikp
    data.push_back(dato(0.09e-2, 0.27e-2, 0.07e-2)); //LHCb:2017joy
    data.push_back(dato(5.9e-3, 3.6e-3, 0.7e-3)); //D0:2013mrm
    data.push_back(dato(-7.6e-3, 5.0e-3, 2.2e-3)); //Belle:2010zqr
    data.push_back(dato(90.e-3, 70.e-3, 20.e-3)); //Belle:2007oni
    data.push_back(dato(30.e-3, 15.e-3, 6.e-3)); //BaBar:2004htr
    data.push_back(dato(18.e-3, 43.e-3, 4.e-3)); //CLEO:2000oig


    pdgaverage.setData(data);
    pdgaverage.setName("ACPBpjpsikp");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    //ACPBpjpsipp
    data.push_back(dato(1.51e-2, 0.50e-2, 0.08e-2)); //LHCb:2024exp
    data.push_back(dato(-4.2e-2, 4.4e-2, 0.9e-2)); //D0:2013mrm
    data.push_back(dato(12.3e-2, 8.5e-2, 0.4e-3)); //BaBar:2004kla
    data.push_back(dato(-2.3e-2, 16.4e-2, 1.5e-2)); //Belle:2002oex


    pdgaverage.setData(data);
    pdgaverage.setName("ACPBpjpsipp");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    //deltaA_Bpjpsipp_Bpjpsikp
    meas.insert(pair<string, dato>("deltaA_Bpjpsipp_Bpjpsikp", dato(1.42e-2, 0.43e-2, 0.08e-2))); //LHCb:2024exp

    //lambda Bsjpsiph from CMS:2020efq
    meas.insert(pair<string, dato>("lambda_Bsjpsiph", dato(0.972, 0.026, 0.008)));

    //CBdjpsikst
    meas.insert(pair<string, dato>("CBdjpsikst", dato(0.025, 0.083, 0.054))); //BaBar:2009byl

    //SBdjpsikst
    meas.insert(pair<string, dato>("SBdjpsikst", dato(0.601, 0.239, 0.087))); //BaBar:2009byl


    //polarization

    //f_0 Bdjpsikst
    data.push_back(dato(0.587, 0.011)); //D0:2008nly
    data.push_back(dato(0.556, 0.009, 0.010)); //BaBar:2007rbr
    data.push_back(dato(0.562, 0.026, 0.018)); //CDF:2004dxr
    data.push_back(dato(0.574, 0.012, 0.009)); //Belle:2005qtf
    data.push_back(dato(0.59, 0.06, 0.01)); //CDF:2000edf
    data.push_back(dato(0.52, 0.07, 0.04)); //CLEO:1997ilq
    data.push_back(dato(0.65, 0.10, 0.04)); //CDF:1995kwt
    data.push_back(dato(0.97, 0.16, 0.15)); //ARGUS:1994rms


    pdgaverage.setData(data);
    pdgaverage.setName("f_0_Bdjpsikst");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    //f_paral Bdjpsikst
    data.push_back(dato(0.230, 0.013, 0.025)); //D0:2008nly
    data.push_back(dato(0.211, 0.010, 0.006)); //BaBar:2007rbr


    pdgaverage.setData(data);
    pdgaverage.setName("f_paral_Bdjpsikst");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    //f_perp Bdjpsikst
    data.push_back(dato(0.233, 0.010, 0.005)); //BaBar:2007rbr
    data.push_back(dato(0.215, 0.032, 0.006)); //CDF:2004dxr
    data.push_back(dato(0.195, 0.012, 0.008)); //Belle:2005qtf


    pdgaverage.setData(data);
    pdgaverage.setName("f_perp_Bdjpsikst");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    //delta_paral Bdjpsikst
    data.push_back(dato(-2.69, 0.08, 0.11)); //D0:2008nly
    data.push_back(dato(-2.93, 0.08, 0.04)); //BaBar:2007rbr




    pdgaverage.setData(data);
    pdgaverage.setName("delta_paral_Bdjpsikst");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    //delta_perp Bdjpsikst
    data.push_back(dato(3.21, 0.06, 0.06)); //D0:2008nly
    data.push_back(dato(2.91, 0.05, 0.03)); //BaBar:2007rbr


    pdgaverage.setData(data);
    pdgaverage.setName("delta_perp_Bdjpsikst");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    //------------------------------------------------

    //Correlated Measurements
    std::vector<dato> CorrData;
    TMatrixDSym CorrStat(2);
    TMatrixDSym CorrSyst(2);


    //Bdjpsik0s : C and S observables from LHCb:2023zcp
    CorrData.push_back(dato(0.010, 0.012));  // C observable
    CorrData.push_back(dato(0.726, 0.014)); // S observable


// Populate the correlation matrix
    CorrStat(0, 0) = 1.;       // Variance for C
    CorrStat(1, 1) = 1.;       // Variance for S
    CorrStat(0, 1) = 0.41;      // Correlation between C and S
    CorrStat(1, 0) = 0.41;      // Symmetric part (same as Corr(0, 1))


// Insert correlated data into corrmeas
    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("CS_Bdjpsik0s_LHCb2023",CorrelatedGaussianObservables(CorrData, CorrStat, CorrSyst)));

    CorrData.clear();


    //Bdjpsik0 : C and S observables from Belle-II:2024lwr
    CorrData.push_back(dato(-0.035, 0.026, 0.029));  // C observable
    CorrData.push_back(dato(0.724, 0.035, 0.009)); // S observable


// Populate the correlation matrix
    CorrStat(0, 0) = 1.;       //
    CorrStat(1, 1) = 1.;
    CorrStat(0, 1) = -0.09;      // Correlation between C and S
    CorrStat(1, 0) = -0.09;      // Symmetric part (same as Corr(0, 1))


    // Populate the correlation matrix
    CorrSyst(0, 0) = 1.;       //
    CorrSyst(1, 1) = 1.;
    CorrSyst(0, 1) = 0.;      // Correlation between C and S
    CorrSyst(1, 0) = 0.;      // Symmetric part (same as Corr(0, 1))


// Insert correlated data into corrmeas
    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("CS_Bdjpsik0s_BelleII2024",CorrelatedGaussianObservables(CorrData, CorrStat, CorrSyst)));

    CorrData.clear();


    //Bdjpsip0 : C and S observables from BaBar:2008kfx
    CorrData.push_back(dato(-0.20, 0.19, 0.03));  // C observable
    CorrData.push_back(dato(-1.23, 0.21, 0.04)); // S observable


// Populate the correlation matrix
    CorrStat(0, 0) = 1.;       //
    CorrStat(1, 1) = 1.;
    CorrStat(0, 1) = 0.197;      // Correlation between C and S
    CorrStat(1, 0) = 0.197;      // Symmetric part (same as Corr(0, 1))


    // Populate the correlation matrix
    CorrSyst(0, 0) = 1.;       //
    CorrSyst(1, 1) = 1.;
    CorrSyst(0, 1) = 0.;      // Correlation between C and S
    CorrSyst(1, 0) = 0.;      // Symmetric part (same as Corr(0, 1))


// Insert correlated data for BaBar:2008kfx into corrmeas
    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("CS_Bdjpsip0_BaBar2008",CorrelatedGaussianObservables(CorrData, CorrStat, CorrSyst)));

    CorrData.clear();

    //C and S for Bdjpsirh from LHCb:2014xpr

    CorrData.push_back(dato(-0.0605, 0.056, 0.0165));  // C observable
    CorrData.push_back(dato(-0.652, 0.125, 0.06)); // S observable


// Populate the correlation matrix
    CorrStat(0, 0) = 1.;       // Variance for C
    CorrStat(1, 1) = 1.;       // Variance for S
    CorrStat(0, 1) = -0.01;      // Correlation between C and S
    CorrStat(1, 0) = -0.01;      // Symmetric part (same as Corr(0, 1))


// Insert correlated data into corrmeas
    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("CS_Bdjpsirh_LHCb2014",CorrelatedGaussianObservables(CorrData, CorrStat, CorrSyst)));

    CorrData.clear();


    //CP violating observables Bsjpsiph from LHCb:2023sim
    TMatrixDSym Corr(6);
    CorrData.push_back(dato(-0.031, 0.018)); //phi_s
    CorrData.push_back(dato(0.990, 0.010)); //lambda
    CorrData.push_back(dato(0.2471, 0.0031)); //f_perp
    CorrData.push_back(dato(0.5175, 0.0035)); //f0
    CorrData.push_back(dato(2.94, 0.07)); //g_perp - g_0
    CorrData.push_back(dato(3.150, 0.062)); //g_par - g_0

    Corr(1,1) = Corr(0,0) = Corr(2,2) = Corr(3,3) = Corr(4,4) = Corr(5,5) = 1.;
    Corr(1,0) = Corr(0,1) = -0.01;
    Corr(2,0) = Corr(0,2) = -0.01;
    Corr(3,0) = Corr(0,3) = -0.01;
    Corr(4,0) = Corr(0,4) = 0.04;
    Corr(5,0) = Corr(0,5) = 0.01;
    Corr(1,2) = Corr(2,1) = 0.00;
    Corr(1,3) = Corr(3,1) = 0.01;
    Corr(1,4) = Corr(4,1) = -0.11;
    Corr(1,5) = Corr(5,1) = -0.03;
    Corr(2,3) = Corr(3,2) = -0.17;
    Corr(2,4) = Corr(4,2) = -0.01;
    Corr(2,5) = Corr(5,2) = -0.04;
    Corr(3,4) = Corr(4,3) = 0.0;
    Corr(3,5) = Corr(5,3) = 0.1;
    Corr(4,5) = Corr(5,4) = 0.2;

    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("phi_lambda_Bsjpsiph_LHCb2023",CorrelatedGaussianObservables(CorrData, Corr)));
    CorrData.clear();

    //phi and lambda for ATLAS:2020lbz solution B
    //TMatrixDSym Corr(5,5);

    CorrData.push_back(dato(-0.081, 0.041, 0.022)); //phi_s
    CorrData.push_back(dato(0.2213, 0.0019, 0.0023)); //f_paral
    CorrData.push_back(dato(0.5131, 0.0013, 0.0038)); //f0
    CorrData.push_back(dato(2.91, 0.11, 0.06)); //g_perp - g_0
    CorrData.push_back(dato(2.94, 0.05, 0.09)); //g_par - g_0

    Corr5(1,1) = Corr5(0,0) = Corr5(2,2) = Corr5(3,3) = Corr5(4,4) = 1.;
    Corr5(1,0) = Corr5(0,1) = -0.003;
    Corr5(2,0) = Corr5(0,2) = -0.004;
    Corr5(3,0) = Corr5(0,3) = 0.004;
    Corr5(4,0) = Corr5(0,4) = 0.007;
    Corr5(1,2) = Corr5(2,1) = -0.341;
    Corr5(1,3) = Corr5(3,1) = 0.133;
    Corr5(1,4) = Corr5(4,1) = 0.522;
    Corr5(2,3) = Corr5(3,2) = -0.034;
    Corr5(2,4) = Corr5(4,2) = -0.103;
    Corr5(3,4) = Corr5(4,3) = 0.254;

    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("phi_Bsjpsiph_ATLAS2020B",CorrelatedGaussianObservables(CorrData, Corr5)));
    CorrData.clear();


    //CP violating observables Bsjpsiph from LHCb:2021wte

    CorrData.push_back(dato(0.00, 0.28, 0.07)); //phi_s
    CorrData.push_back(dato(0.885, 0.014, 0.031)); //lambda
    CorrData.push_back(dato(0.234, 0.034, 0.008)); //f_perp
    CorrData.push_back(dato(0.530, 0.029, 0.013)); //f0
    CorrData.push_back(dato(2.415, 0.425, 0.100)); //g_perp - g_0
    CorrData.push_back(dato(3.115, 0.075, 0.060)); //g_par - g_0

    Corr(1,1) = Corr(0,0) = Corr(2,2) = Corr(3,3) = Corr(4,4) = Corr(5,5) = 1.;
    Corr(1,0) = Corr(0,1) = 0.15;
    Corr(2,0) = Corr(0,2) = -0.06;
    Corr(3,0) = Corr(0,3) = 0.07;
    Corr(4,0) = Corr(0,4) = 0.08;
    Corr(5,0) = Corr(0,5) = 0.13;
    Corr(1,2) = Corr(2,1) = -0.14;
    Corr(1,3) = Corr(3,1) = 0.24;
    Corr(1,4) = Corr(4,1) = -0.11;
    Corr(1,5) = Corr(5,1) = -0.06;
    Corr(2,3) = Corr(3,2) = -0.66;
    Corr(2,4) = Corr(4,2) = 0.10;
    Corr(2,5) = Corr(5,2) = -0.06;
    Corr(3,4) = Corr(4,3) = -0.17;
    Corr(3,5) = Corr(5,3) = 0.08;
    Corr(4,5) = Corr(5,4) = -0.03;

    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("phi_lambda_Bsjpsiph_LHCb2021",CorrelatedGaussianObservables(CorrData, Corr)));
    CorrData.clear();

    //CP violating observables Bsjpsiph from CMS:2020efq
    //TMatrixDSym Corr5(5,5);

    CorrData.push_back(dato(-0.021, 0.044, 0.010)); //phi_s
    CorrData.push_back(dato(0.2393, 0.0050, 0.0037)); //f_perp
    CorrData.push_back(dato(0.5289, 0.0038, 0.0041)); //f0
    CorrData.push_back(dato(2.78, 0.15, 0.06)); //g_perp - g_0
    CorrData.push_back(dato(3.19, 0.12, 0.04)); //g_par - g_0

    Corr5(1,1) = Corr5(0,0) = Corr5(2,2) = Corr5(3,3) = Corr5(4,4) = 1.;
    Corr5(1,0) = Corr5(0,1) = -0.01;
    Corr5(2,0) = Corr5(0,2) = 0.01;
    Corr5(3,0) = Corr5(0,3) = -0.08;
    Corr5(4,0) = Corr5(0,4) = -0.01;
    Corr5(1,2) = Corr5(2,1) = -0.56;
    Corr5(1,3) = Corr5(3,1) = 0.01;
    Corr5(1,4) = Corr5(4,1) = -0.03;
    Corr5(2,3) = Corr5(3,2) = 0.01;
    Corr5(2,4) = Corr5(4,2) = -0.03;
    Corr5(3,4) = Corr5(4,3) = 0.26;


    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("phi_Bsjpsiph_CMS2020",CorrelatedGaussianObservables(CorrData, Corr5)));
    CorrData.clear();


    //phi and lambda for Bsjpsiph from LHCb:2019sgv
    CorrData.push_back(dato(0.002, 0.044, 0.012));  // phi
    CorrData.push_back(dato(0.949, 0.036, 0.019)); // lambda


// Populate the correlation matrix
    CorrStat(0, 0) = 1.;       //
    CorrStat(1, 1) = 1.;
    CorrStat(0, 1) = 0.026;      // Correlation between C and S
    CorrStat(1, 0) = 0.026;      // Symmetric part (same as Corr(0, 1))


    // Populate the correlation matrix
    CorrSyst(0, 0) = 1.;       //
    CorrSyst(1, 1) = 1.;
    CorrSyst(0, 1) = 0.;      // Correlation between C and S
    CorrSyst(1, 0) = 0.;      // Symmetric part (same as Corr(0, 1))

    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("phi_lambda_Bsjpsiph_LHCb2019",CorrelatedGaussianObservables(CorrData, CorrStat, CorrSyst)));

    CorrData.clear();

    //phi and lmbda for Bsjpsiph from LHCb:2017hbp
    CorrData.push_back(dato(-0.025, 0.044, 0.008));  // phi
    CorrData.push_back(dato(0.978, 0.013, 0.003)); // lambda


// Populate the correlation matrix
    CorrStat(0, 0) = 1.;       //
    CorrStat(1, 1) = 1.;
    CorrStat(0, 1) = -0.04;      // Correlation between C and S
    CorrStat(1, 0) = -0.04;      // Symmetric part (same as Corr(0, 1))


    // Populate the correlation matrix
    CorrSyst(0, 0) = 1.;       //
    CorrSyst(1, 1) = 1.;
    CorrSyst(0, 1) = 0.;      // Correlation between C and S
    CorrSyst(1, 0) = 0.;      // Symmetric part (same as Corr(0, 1))

    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("phi_lambda_Bsjpsiph_LHCb2017",CorrelatedGaussianObservables(CorrData, CorrStat, CorrSyst)));

    CorrData.clear();

    //CP violating observables Bsjpsiph from ATLAS:2016pno

    CorrData.push_back(dato(-0.110, 0.082, 0.042)); //phi_s
    CorrData.push_back(dato(0.230, 0.005, 0.006)); //f_paral
    CorrData.push_back(dato(0.520, 0.004, 0.007)); //f0
    CorrData.push_back(dato(4.50, 0.45, 0.30)); //g_perp - g_0
    CorrData.push_back(dato(3.15, 0.10, 0.05)); //g_par - g_0

    Corr5(1,1) = Corr5(0,0) = Corr5(2,2) = Corr5(3,3) = Corr5(4,4) = 1.;
    Corr5(1,0) = Corr5(0,1) = 0.030;
    Corr5(2,0) = Corr5(0,2) = 0.029;
    Corr5(3,0) = Corr5(0,3) = 0.035;
    Corr5(4,0) = Corr5(0,4) = 0.067;
    Corr5(1,2) = Corr5(2,1) = -0.330;
    Corr5(1,3) = Corr5(3,1) = 0.105;
    Corr5(1,4) = Corr5(4,1) = -0.03;
    Corr5(2,3) = Corr5(3,2) = 0.007;
    Corr5(2,4) = Corr5(4,2) = -0.011;
    Corr5(3,4) = Corr5(4,3) = 0.158;


    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("phi_Bsjpsiph_ATLAS2016",CorrelatedGaussianObservables(CorrData, Corr5)));
    CorrData.clear();

    //CP violating observables Bsjpsiph from ATLAS:2014nmm
    //TMatrixDSym Corr4(4);
    CorrData.push_back(dato(0.12, 0.25, 0.05)); //phi_s
    CorrData.push_back(dato(0.220, 0.008, 0.009)); //f_paral
    CorrData.push_back(dato(0.529, 0.006, 0.012)); //f0
    CorrData.push_back(dato(3.89, 0.47, 0.11)); //g_perp - g_0
    CorrData.push_back(dato(3.135, 0.095, 0.09)); //g_paral - g_0


    Corr5(1,1) = Corr5(0,0) = Corr5(2,2) = Corr5(3,3) = Corr(5,5) = 1.;
    Corr5(1,0) = Corr5(0,1) = 0.010;
    Corr5(2,0) = Corr5(0,2) = 0.002;
    Corr5(3,0) = Corr5(0,3) = -0.043;
    Corr5(4,0) = Corr5(0,4) = 0.021;
    Corr5(1,2) = Corr5(2,1) = -0.316;
    Corr5(1,3) = Corr5(3,1) = 0.005;
    Corr5(1,4) = Corr5(4,1) = 0.008;
    Corr5(2,3) = Corr5(3,2) = -0.016;
    Corr5(2,4) = Corr5(4,2) = -0.003;
    Corr5(3,4) = Corr5(4,3) = 0.038;


    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("phi_Bsjpsiph_ATLAS2014",CorrelatedGaussianObservables(CorrData, Corr5)));
    CorrData.clear();


    //polarizations in Bdjpsikst from LHCb:2013vga
    TMatrixDSym Corr4(4);
    CorrData.push_back(dato(0.227, 0.004, 0.011)); //f_paral
    CorrData.push_back(dato(0.201, 0.004, 0.008)); //f_perp
    CorrData.push_back(dato(-2.94, 0.02, 0.03)); //delta_par
    CorrData.push_back(dato(2.94, 0.02, 0.02)); //delta_perp

    Corr4(1,1) = Corr4(0,0) = Corr4(2,2) = Corr4(3,3) = 1.;
    Corr4(1,0) = Corr4(0,1) = -0.70;
    Corr4(2,0) = Corr4(0,2) = 0.12;
    Corr4(3,0) = Corr4(0,3) = 0.04;
    Corr4(1,2) = Corr4(2,1) = -0.14;
    Corr4(1,3) = Corr4(3,1) = -0.01;
    Corr4(2,3) = Corr4(3,2) = 0.64;


    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("polarization_Bdjpsikst_LHCb2013",CorrelatedGaussianObservables(CorrData, Corr4)));
    CorrData.clear();


    //Cp polarization obs for Bsjpsikst from LHCb:2015esn

    TMatrixDSym Corr7(7);
    CorrData.push_back(dato(-0.048, 0.057, 0.020)); //A0 CP
    CorrData.push_back(dato(0.171, 0.152, 0.028)); //A_paral CP
    CorrData.push_back(dato(-0.049, 0.096, 0.025)); //A_perp CP
    CorrData.push_back(dato(0.497, 0.025, 0.025)); // f0
    CorrData.push_back(dato(0.179, 0.027, 0.013)); // f_paral
    CorrData.push_back(dato(-2.70, 0.16, 0.013)); //delta_paral
    CorrData.push_back(dato(0.0105, 0.11, 0.0165)); //delta_perp

    Corr7(1,1) = Corr7(0,0) = Corr7(2,2) = Corr7(3,3) = Corr7(4,4) = Corr7(5,5) = Corr7(6,6) = 1.;
    Corr7(1,0) = Corr7(0,1) = -0.11;
    Corr7(2,0) = Corr7(0,2) = -0.17;
    Corr7(3,0) = Corr7(0,3) = 0.06;
    Corr7(4,0) = Corr7(0,4) = -0.05;
    Corr7(5,0) = Corr7(0,5) = 0.03;
    Corr7(6,0) = Corr7(0,6) = 0.02;
    Corr7(1,2) = Corr7(2,1) = -0.49;
    Corr7(1,3) = Corr7(3,1) = -0.04;
    Corr7(1,4) = Corr7(4,1) = -0.07;
    Corr7(1,5) = Corr7(5,1) = 0.09;
    Corr7(1,6) = Corr7(6,1) = 0.06;
    Corr7(2,3) = Corr7(3,2) = 0.01;
    Corr7(2,4) = Corr7(4,2) = -0.06;
    Corr7(2,5) = Corr7(5,2) = -0.09;
    Corr7(2,6) = Corr7(6,2) = -0.03;
    Corr7(3,4) = Corr7(4,3) = -0.34;
    Corr7(3,5) = Corr7(5,3) = 0.04;
    Corr7(3,6) = Corr7(6,3) = 0.05;
    Corr7(4,5) = Corr7(5,4) = -0.03;
    Corr7(4,6) = Corr7(6,4) = -0.04;
    Corr7(5,6) = Corr7(6,5) = 0.62;

    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("polarization_Bsjpsikst_LHCb2015",CorrelatedGaussianObservables(CorrData, Corr)));
    CorrData.clear();

    //data for the pull in newmeas

    //CBdjpsik0
    data.push_back(dato(0.010, 0.012));
    data.push_back(dato(-0.035, 0.026, 0.029));
    data.push_back(dato(0.005, 0.021, 0.037));

    pdgaverage.setData(data);
    pdgaverage.setName("CBdjpsik0s");
    pdgaverage.CalculateAverage();

    newmeas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    //SBdjpsik0

    data.push_back(dato(0.726, 0.014));
    data.push_back(dato(0.724, 0.035, 0.009));
    data.push_back(dato(0.670, 0.029, 0.013));  ////Belle:2012paq
    data.push_back(dato(0.57, 0.58, 0.06)); //Belle:2012teq
    data.push_back(dato(0.775, 0.425)); //CDF:1999ijp
    data.push_back(dato(0.73, 0.93, 0.16)); //ALEPH:2000jem
    data.push_back(dato(3.1, 1.9, 0.5)); //OPAL:1998mtm

    pdgaverage.setData(data);
    pdgaverage.setName("SBdjpsik0s");
    pdgaverage.CalculateAverage();

    newmeas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    //CBdjpsip0

    data.push_back(dato(-0.20, 0.19, 0.03));
    data.push_back(dato(0.155, 0.14, 0.035));  ////Belle:2018nxw
    data.push_back(dato(0.13, 0.12, 0.03)); //Belle-II:2024hqw

    pdgaverage.setData(data);
    pdgaverage.setName("CBdjpsip0");
    pdgaverage.CalculateAverage();

    newmeas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    //SBdjpsip0

    data.push_back(dato(-1.23, 0.21, 0.04));
    data.push_back(dato(-0.59, 0.19, 0.03));  ////Belle:2018nxw
    data.push_back(dato(-0.88, 0.17, 0.03)); //Belle-II:2024hqw

    pdgaverage.setData(data);
    pdgaverage.setName("SBdjpsip0");
    pdgaverage.CalculateAverage();

    newmeas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    //phi BsJpsiph
    data.push_back(dato(-0.031, 0.018)); //LHCb:2023sim
    data.push_back(dato(-0.021, 0.044, 0.010)); //CMS:2020efq
    data.push_back(dato(-0.081, 0.041, 0.022)); //ATLAS:2020lbz
    data.push_back(dato(0.00, 0.28, 0.07)); //LHCb:2021wte
    data.push_back(dato(-0.110, 0.082, 0.042)); //ATLAS:2016pno
    data.push_back(dato(0.002, 0.044, 0.012));
    data.push_back(dato(-0.025, 0.044, 0.008));
    data.push_back(dato(0.12, 0.25, 0.05)); //ATLAS:2014nmm


    pdgaverage.setData(data);
    pdgaverage.setName("phi_Bsjpsiph");
    pdgaverage.CalculateAverage();

    newmeas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();


    //lambda BsJpsiph
    data.push_back(dato(0.990, 0.010)); //LHCb:2023sim
    data.push_back(dato(0.885, 0.014, 0.031)); //LHCb:2021wte
    data.push_back(dato(0.972, 0.026, 0.008)); //CMS:2020efq
    data.push_back(dato(0.949, 0.036, 0.019));
    data.push_back(dato(0.978, 0.013, 0.003));



    pdgaverage.setData(data);
    pdgaverage.setName("lambda_Bsjpsiph");
    pdgaverage.CalculateAverage();

    newmeas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();


    //f_perp
    data.push_back(dato(0.2471, 0.0031));
    data.push_back(dato(0.234, 0.034, 0.008));
    data.push_back(dato(0.2393, 0.0050, 0.0037));



    pdgaverage.setData(data);
    pdgaverage.setName("f_perp_Bsjpsiph");
    pdgaverage.CalculateAverage();

    newmeas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    //f_paral
    data.push_back(dato(0.2213, 0.0019, 0.0023));
    data.push_back(dato(0.230, 0.005, 0.006));
    data.push_back(dato(0.220, 0.008, 0.009));



    pdgaverage.setData(data);
    pdgaverage.setName("f_paral_Bsjpsiph");
    pdgaverage.CalculateAverage();

    newmeas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    //f_0
    data.push_back(dato(0.5175, 0.0035)); //LHCb:2023sim
    data.push_back(dato(0.5131, 0.0013, 0.0038));
    data.push_back(dato(0.530, 0.029, 0.013));
    data.push_back(dato(0.5289, 0.0038, 0.0041));
    data.push_back(dato(0.520, 0.004, 0.007));
    data.push_back(dato(0.529, 0.006, 0.012));


    pdgaverage.setData(data);
    pdgaverage.setName("f_0_Bsjpsiph");
    pdgaverage.CalculateAverage();

    newmeas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    //delta_perp
    data.push_back(dato(2.94, 0.07)); //LHCb:2023sim
    data.push_back(dato(2.91, 0.11, 0.06));
    data.push_back(dato(2.415, 0.425, 0.100));
    data.push_back(dato(2.78, 0.15, 0.06));
    data.push_back(dato(4.50, 0.45, 0.30));
    data.push_back(dato(3.89, 0.47, 0.11));


    pdgaverage.setData(data);
    pdgaverage.setName("delta_perp_Bsjpsiph");
    pdgaverage.CalculateAverage();

    newmeas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    //delta_paral
    data.push_back(dato(3.150, 0.062)); //LHCb:2023sim
    data.push_back(dato(2.94, 0.05, 0.09));
    data.push_back(dato(3.115, 0.075, 0.060));
    data.push_back(dato(3.19, 0.12, 0.04));
    data.push_back(dato(3.15, 0.10, 0.05));
    data.push_back(dato(3.89, 0.47, 0.11));
    data.push_back(dato(3.135, 0.095, 0.09));

    pdgaverage.setData(data);
    pdgaverage.setName("delta_paral_Bsjpsiph");
    pdgaverage.CalculateAverage();

    newmeas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    //fparal Bdjpsikst
    data.push_back(dato(0.227, 0.004, 0.011));
    data.push_back(dato(0.230, 0.013, 0.025)); //D0:2008nly
    data.push_back(dato(0.211, 0.010, 0.006)); //BaBar:2007rbr


    pdgaverage.setData(data);
    pdgaverage.setName("f_paral_Bdjpsikst");
    pdgaverage.CalculateAverage();

    newmeas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();


    //fperp Bdjpsikst
    data.push_back(dato(0.233, 0.010, 0.005)); //BaBar:2007rbr
    data.push_back(dato(0.215, 0.032, 0.006)); //CDF:2004dxr
    data.push_back(dato(0.195, 0.012, 0.008)); //Belle:2005qtf
    data.push_back(dato(0.201, 0.004, 0.008));


    pdgaverage.setData(data);
    pdgaverage.setName("f_perp_Bdjpsikst");
    pdgaverage.CalculateAverage();

    newmeas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    //delta_paral Bdjpsikst
    data.push_back(dato(-2.69, 0.08, 0.11)); //D0:2008nly
    data.push_back(dato(-2.93, 0.08, 0.04)); //BaBar:2007rbr
    data.push_back(dato(-2.94, 0.02, 0.03));


    pdgaverage.setData(data);
    pdgaverage.setName("delta_paral_Bdjpsikst");
    pdgaverage.CalculateAverage();

    newmeas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    //delta_perp Bdjpsikst
    data.push_back(dato(3.21, 0.06, 0.06)); //D0:2008nly
    data.push_back(dato(2.91, 0.05, 0.03)); //BaBar:2007rbr
    data.push_back(dato(2.94, 0.02, 0.02));


    pdgaverage.setData(data);
    pdgaverage.setName("delta_perp_Bdjpsikst");
    pdgaverage.CalculateAverage();

    newmeas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    std::cout << "All meas inserted" << endl;


   // Create histograms for uncorrelated observables (from meas)
   histos.createH1D("R_Bpjpsipp_Bpjpsikp", 500, 0.0, 0.0);
   histos.createH1D("R_Bdjpsiom_Bdjpsirh", 500, 0.0, 0.0);
   histos.createH1D("R_Bdjpsikst_Bdjpsik0", 500, 0.0, 0.0);
   histos.createH1D("deltaA_Bpjpsipp_Bpjpsikp", 500, 0.0, 0.0);
    for (const auto& channel : channelNamesSU3) {
      //create histos for mod and phase as well as real and imaginary parts for the effective parameters
      if (std::find(vectorMesonChannels.begin(), vectorMesonChannels.end(), channel) != vectorMesonChannels.end()) {

        // Vector meson channels: store histograms for all polarization components
        histos.createH1D("mod_A_0_" + channel, 500, 0.0, 0.0);
        histos.createH1D("phase_A_0_" + channel, 500, 0.0, 0.0);
        histos.createH2D("mod_A_0_" + channel, "phase_A_0_" + channel, 500, 0.0, 0.0, 500, 0.0, 0.0);
        histos.createH1D("A_0_re_" + channel, 500, 0.0, 0.0);
        histos.createH1D("A_0_im_" + channel, 500, 0.0, 0.0);

        histos.createH1D("mod_B_0_" + channel, 500, 0.0, 0.0);
        histos.createH1D("phase_B_0_" + channel, 500, 0.0, 0.0);
        histos.createH2D("mod_B_0_" + channel, "phase_B_0_" + channel, 500, 0.0, 0.0, 500, 0.0, 0.0);
        histos.createH1D("B_0_re_" + channel, 500, 0.0, 0.0);
        histos.createH1D("B_0_im_" + channel, 500, 0.0, 0.0);

        // Perpendicular polarization
        histos.createH1D("mod_A_perp_" + channel, 500, 0.0, 0.0);
        histos.createH1D("phase_A_perp_" + channel, 500, 0.0, 0.0);
        histos.createH2D("mod_A_perp_" + channel, "phase_A_perp_" + channel, 500, 0.0, 0.0, 500, 0.0, 0.0);
        histos.createH1D("A_perp_re_" + channel, 500, 0.0, 0.0);
        histos.createH1D("A_perp_im_" + channel, 500, 0.0, 0.0);

        histos.createH1D("mod_B_perp_" + channel, 500, 0.0, 0.0);
        histos.createH1D("phase_B_perp_" + channel, 500, 0.0, 0.0);
        histos.createH2D("mod_B_perp_" + channel, "phase_B_perp_" + channel, 500, 0.0, 0.0, 500, 0.0, 0.0);
        histos.createH1D("B_perp_re_" + channel, 500, 0.0, 0.0);
        histos.createH1D("B_perp_im_" + channel, 500, 0.0, 0.0);

        // Parallel polarization
        histos.createH1D("mod_A_paral_" + channel, 500, 0.0, 0.0);
        histos.createH1D("phase_A_paral_" + channel, 500, 0.0, 0.0);
        histos.createH2D("mod_A_paral_" + channel, "phase_A_paral_" + channel, 500, 0.0, 0.0, 500, 0.0, 0.0);
        histos.createH1D("A_paral_re_" + channel, 500, 0.0, 0.0);
        histos.createH1D("A_paral_im_" + channel, 500, 0.0, 0.0);

        histos.createH1D("mod_B_paral_" + channel, 500, 0.0, 0.0);
        histos.createH1D("phase_B_paral_" + channel, 500, 0.0, 0.0);
        histos.createH2D("mod_B_paral_" + channel, "phase_B_paral_" + channel, 500, 0.0, 0.0, 500, 0.0, 0.0);
        histos.createH1D("B_paral_re_" + channel, 500, 0.0, 0.0);
        histos.createH1D("B_paral_im_" + channel, 500, 0.0, 0.0);

    } else {
        // Pseudoscalar channels: only A and B
        histos.createH1D("mod_A_" + channel, 500, 0.0, 0.0);
        histos.createH1D("phase_A_" + channel, 500, 0.0, 0.0);
        histos.createH2D("mod_A_" + channel, "phase_A_" + channel, 500, 0.0, 0.0, 500, 0.0, 0.0);
        histos.createH1D("A_re_" + channel, 500, 0.0, 0.0);
        histos.createH1D("A_im_" + channel, 500, 0.0, 0.0);

        histos.createH1D("mod_B_" + channel, 500, 0.0, 0.0);
        histos.createH1D("phase_B_" + channel, 500, 0.0, 0.0);
        histos.createH2D("mod_B_" + channel, "phase_B_" + channel, 500, 0.0, 0.0, 500, 0.0, 0.0);
        histos.createH1D("B_re_" + channel, 500, 0.0, 0.0);
        histos.createH1D("B_im_" + channel, 500, 0.0, 0.0);
    }
         // Branching ratio (BR)
        if (meas.find("BR" + channel) != meas.end()) {
            if (histos.h1d.find("BR_" + channel) == histos.h1d.end()) {  // Check if it already exists
                histos.createH1D("BR_" + channel, 500, 0.0, 0.0);

            }
        }

        // ACP, C, and S observables
        if (meas.find("ACP" + channel) != meas.end()) {
            if (histos.h1d.find("ACP_" + channel) == histos.h1d.end()) {
                histos.createH1D("ACP_" + channel, 500, 0., 0.);

            }
        }
        if (meas.find("C" + channel) != meas.end()) {
            if (histos.h1d.find("C_" + channel) == histos.h1d.end()) {
                histos.createH1D("C_" + channel, 500, 0.0, 0.0);

            }
        }
        if (meas.find("S" + channel) != meas.end()) {
            if (histos.h1d.find("S_" + channel) == histos.h1d.end()) {
                histos.createH1D("S_" + channel, 500, 0.0, 0.0);

            }
        }

    }

    // Create histograms for correlated observables (from corrmeas)
    for (auto& pair : corrmeas) {
        const auto& corr_key = pair.first;
        auto& corrObservable = pair.second;

        // Parse the channels from the key
        auto channel_and_experiment = extractChannelFromCorrKey(corr_key);
        auto& channels = channel_and_experiment.first;

        // Handle C and S observables (e.g., "CS_" prefix)
        if (corr_key.find("CS_") != std::string::npos) {
            if (histos.h1d.find("C_" + channels[0]) == histos.h1d.end()) {
                histos.createH1D("C_" + channels[0], 500, 0.0, 0.0);

            }
            if (histos.h1d.find("S_" + channels[0]) == histos.h1d.end()) {
                histos.createH1D("S_" + channels[0], 500, 0.0, 0.0);

            }
        }

        // Handle ACP observables (e.g., "ACP_" prefix)
        if (corr_key.find("ACP_") != std::string::npos) {
            for (size_t i = 0; i < channels.size(); ++i) {
                if (histos.h1d.find("ACP_" + channels[i]) == histos.h1d.end()) {
                    histos.createH1D("ACP_" + channels[i], 500, 0.0, 0.0);

                }
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

  //method to define the parameters needed to calculate each decay amplitude, uses BCModel method AddParameter
void goldenmodes::DefineParameters(const string& channel) {
   double limit1 = 0.;


  if (channel == "Bdjpsik0") {
        std::vector<std::string> params = {
	  "A_Bdjpsik0_re",
	  "A_Bdjpsik0_im",
	  "B_Bdjpsik0_re",
	  "B_Bdjpsik0_im"
        };
        channelParameters[channel] = params;
        for (const auto& param : params) {
            AddParameter(param, -50., 50.);
        }
	GetParameter("A_Bdjpsik0_re").SetLowerLimit(0.);
	GetParameter("A_Bdjpsik0_im").SetUpperLimit(0.);
	GetParameter("A_Bdjpsik0_im").SetLowerLimit(0.);
  GetParameter("A_Bdjpsik0_re").SetUpperLimit(100.);

} else if (channel == "Bdjpsip0") {
        std::vector<std::string> params = {
      "A_Bdjpsip0_re",
	    "A_Bdjpsip0_im",
      "B_Bdjpsip0_re",
	    "B_Bdjpsip0_im"
    };
    channelParameters[channel] = params;
    for (const auto& param : params) {
      AddParameter(param, -100., 100.);
    }
	GetParameter("A_Bdjpsip0_re").SetLowerLimit(0.);
	//GetParameter("A_Bdjpsip0_re").SetUpperLimit(100.);
  GetParameter("A_Bdjpsip0_im").SetUpperLimit(0.);
	GetParameter("A_Bdjpsip0_im").SetLowerLimit(0.);
} else if (channel == "Bdjpsiom") {
    std::vector<std::string> params = {
      "A_Bdjpsiom_0_re",
	    "A_Bdjpsiom_0_im",
      "B_Bdjpsiom_0_re",
	    "B_Bdjpsiom_0_im",
      "A_Bdjpsiom_perp_re",
	    "A_Bdjpsiom_perp_im",
      "B_Bdjpsiom_perp_re",
	    "B_Bdjpsiom_perp_im",
      "A_Bdjpsiom_paral_re",
	    "A_Bdjpsiom_paral_im",
      "B_Bdjpsiom_paral_re",
	    "B_Bdjpsiom_paral_im"
        };
        channelParameters[channel] = params;
        for (const auto& param : params) {
            AddParameter(param, -50., 50.);
        }
	GetParameter("A_Bdjpsiom_0_re").SetLowerLimit(0.);
	GetParameter("A_Bdjpsiom_0_im").SetUpperLimit(0.);
	GetParameter("A_Bdjpsiom_0_im").SetLowerLimit(0.);
  GetParameter("A_Bdjpsiom_paral_re").SetLowerLimit(0.);
	GetParameter("A_Bdjpsiom_paral_im").SetUpperLimit(0.);
	GetParameter("A_Bdjpsiom_paral_im").SetLowerLimit(0.);
  GetParameter("A_Bdjpsiom_perp_re").SetLowerLimit(0.);
	GetParameter("A_Bdjpsiom_perp_im").SetUpperLimit(0.);
	GetParameter("A_Bdjpsiom_perp_im").SetLowerLimit(0.);
} else if (channel == "Bpjpsikp") {
        std::vector<std::string> params = {
      "A_Bpjpsikp_re",
	    "A_Bpjpsikp_im",
      "B_Bpjpsikp_re",
	    "B_Bpjpsikp_im"
    };
    channelParameters[channel] = params;
    for (const auto& param : params) {
      AddParameter(param, -100., 100.);
    }
    GetParameter("A_Bpjpsikp_im").SetUpperLimit(0.);
	  GetParameter("A_Bpjpsikp_re").SetLowerLimit(0.);
    GetParameter("A_Bpjpsikp_im").SetLowerLimit(0.);
  } else if (channel == "Bpjpsipp") {
    std::vector<std::string> params = {
      "A_Bpjpsipp_re",
	    "A_Bpjpsipp_im",
      "B_Bpjpsipp_re",
	    "B_Bpjpsipp_im"
    };
    channelParameters[channel] = params;
    for (const auto& param : params) {
        AddParameter(param, -50., 50.);
    }
    GetParameter("A_Bpjpsipp_im").SetUpperLimit(0.);
    GetParameter("A_Bpjpsipp_im").SetLowerLimit(0.);
    GetParameter("A_Bpjpsipp_re").SetLowerLimit(0.);
  } else if (channel == "Bsjpsiph") {
      std::vector<std::string> params = {
        "A_Bsjpsiph_0_re",
  	    "A_Bsjpsiph_0_im",
        "B_Bsjpsiph_0_re",
  	    "B_Bsjpsiph_0_im",
        "A_Bsjpsiph_perp_re",
  	    "A_Bsjpsiph_perp_im",
        "B_Bsjpsiph_perp_re",
  	    "B_Bsjpsiph_perp_im",
        "A_Bsjpsiph_paral_re",
  	    "A_Bsjpsiph_paral_im",
        "B_Bsjpsiph_paral_re",
  	    "B_Bsjpsiph_paral_im"
      };
      channelParameters[channel] = params;
      for (const auto& param : params) {
          AddParameter(param, -50., 50.);
      }
      GetParameter("A_Bsjpsiph_0_im").SetUpperLimit(0.);
      GetParameter("A_Bsjpsiph_0_im").SetLowerLimit(0.);
      GetParameter("A_Bsjpsiph_0_re").SetLowerLimit(0.);
      GetParameter("A_Bsjpsiph_paral_im").SetUpperLimit(0.);
      GetParameter("A_Bsjpsiph_paral_im").SetLowerLimit(0.);
      GetParameter("A_Bsjpsiph_paral_re").SetLowerLimit(0.);
      GetParameter("A_Bsjpsiph_perp_im").SetUpperLimit(0.);
      GetParameter("A_Bsjpsiph_perp_im").SetLowerLimit(0.);
      GetParameter("A_Bsjpsiph_perp_re").SetLowerLimit(0.);

    } else if (channel == "Bsjpsik0") {
      std::vector<std::string> params = {
        "A_Bsjpsik0_re",
        "A_Bsjpsik0_im",
        "B_Bsjpsik0_re",
        "B_Bsjpsik0_im"
      };
      channelParameters[channel] = params;
      for (const auto& param : params) {
        AddParameter(param, -50., 50.);
      }
        GetParameter("A_Bsjpsik0_im").SetUpperLimit(0.);
        GetParameter("A_Bsjpsik0_im").SetLowerLimit(0.);
        GetParameter("A_Bsjpsik0_re").SetLowerLimit(0.);
      } else if (channel == "Bdjpsikst") {
        std::vector<std::string> params = {
          "A_Bdjpsikst_0_re",
    	    "A_Bdjpsikst_0_im",
          "B_Bdjpsikst_0_re",
    	    "B_Bdjpsikst_0_im",
          "A_Bdjpsikst_perp_re",
    	    "A_Bdjpsikst_perp_im",
          "B_Bdjpsikst_perp_re",
    	    "B_Bdjpsikst_perp_im",
          "A_Bdjpsikst_paral_re",
    	    "A_Bdjpsikst_paral_im",
          "B_Bdjpsikst_paral_re",
    	    "B_Bdjpsikst_paral_im"
            };
            channelParameters[channel] = params;
            for (const auto& param : params) {
                AddParameter(param, -50., 50.);
            }
    	GetParameter("A_Bdjpsikst_0_re").SetLowerLimit(0.);
    	GetParameter("A_Bdjpsikst_0_im").SetUpperLimit(0.);
    	GetParameter("A_Bdjpsikst_0_im").SetLowerLimit(0.);
      GetParameter("A_Bdjpsikst_paral_re").SetLowerLimit(0.);
    	GetParameter("A_Bdjpsikst_paral_im").SetUpperLimit(0.);
    	GetParameter("A_Bdjpsikst_paral_im").SetLowerLimit(0.);
      GetParameter("A_Bdjpsikst_perp_re").SetLowerLimit(0.);
    	GetParameter("A_Bdjpsikst_perp_im").SetUpperLimit(0.);
    	GetParameter("A_Bdjpsikst_perp_im").SetLowerLimit(0.);
    } else if (channel == "Bdjpsirh") {
      std::vector<std::string> params = {
        "A_Bdjpsirh_0_re",
        "A_Bdjpsirh_0_im",
        "B_Bdjpsirh_0_re",
        "B_Bdjpsirh_0_im",
        "A_Bdjpsirh_perp_re",
        "A_Bdjpsirh_perp_im",
        "B_Bdjpsirh_perp_re",
        "B_Bdjpsirh_perp_im",
        "A_Bdjpsirh_paral_re",
        "A_Bdjpsirh_paral_im",
        "B_Bdjpsirh_paral_re",
        "B_Bdjpsirh_paral_im"
          };
          channelParameters[channel] = params;
          for (const auto& param : params) {
              AddParameter(param, -50., 50.);
          }
    GetParameter("A_Bdjpsirh_0_re").SetLowerLimit(0.);
    GetParameter("A_Bdjpsirh_0_im").SetUpperLimit(0.);
    GetParameter("A_Bdjpsirh_0_im").SetLowerLimit(0.);
    GetParameter("A_Bdjpsirh_paral_re").SetLowerLimit(0.);
    GetParameter("A_Bdjpsirh_paral_im").SetUpperLimit(0.);
    GetParameter("A_Bdjpsirh_paral_im").SetLowerLimit(0.);
    GetParameter("A_Bdjpsirh_perp_re").SetLowerLimit(0.);
    GetParameter("A_Bdjpsirh_perp_im").SetUpperLimit(0.);
    GetParameter("A_Bdjpsirh_perp_im").SetLowerLimit(0.);
    } else if (channel == "Bsjpsikst") {
      std::vector<std::string> params = {
        "A_Bsjpsikst_0_re",
        "A_Bsjpsikst_0_im",
        "B_Bsjpsikst_0_re",
        "B_Bsjpsikst_0_im",
        "A_Bsjpsikst_perp_re",
        "A_Bsjpsikst_perp_im",
        "B_Bsjpsikst_perp_re",
        "B_Bsjpsikst_perp_im",
        "A_Bsjpsikst_paral_re",
        "A_Bsjpsikst_paral_im",
        "B_Bsjpsikst_paral_re",
        "B_Bsjpsikst_paral_im"
          };
          channelParameters[channel] = params;
          for (const auto& param : params) {
              AddParameter(param, -50., 50.);
          }
    GetParameter("A_Bsjpsikst_0_re").SetLowerLimit(0.);
    GetParameter("A_Bsjpsikst_0_im").SetUpperLimit(0.);
    GetParameter("A_Bsjpsikst_0_im").SetLowerLimit(0.);
    GetParameter("A_Bsjpsikst_paral_re").SetLowerLimit(0.);
    GetParameter("A_Bsjpsikst_paral_im").SetUpperLimit(0.);
    GetParameter("A_Bsjpsikst_paral_im").SetLowerLimit(0.);
    GetParameter("A_Bsjpsikst_perp_re").SetLowerLimit(0.);
    GetParameter("A_Bsjpsikst_perp_im").SetUpperLimit(0.);
    GetParameter("A_Bsjpsikst_perp_im").SetLowerLimit(0.);
  } else {
        // Handle the case where the channel is not recognized
        std::cerr << "Error: Unrecognized channel \"" << channel << "\" in DefineParameters." << std::endl;
        throw std::runtime_error("Unrecognized channel: " + channel);
	}

  std::cout << "Number of parameters: " << GetNParameters() << std::endl;
  SetPriorConstantAll();

}


std::map<std::string, double> parameterValues;
 //function that given the channels and the string inside the map channelParameters adds every parameter name to a map <string, double> and returns said map. now parameterValues[channelParameters] is a complex double.
std::map<std::string, double> goldenmodes::DeclareParameters() {

    // Ensure channelParameters is populated
    for (const auto& channel : channelNamesSU3) {
        // Call DefineParameters once per channel
        if (channelParameters.find(channel) == channelParameters.end()) {
            goldenmodes::DefineParameters(channel);
        }

        // Add parameters to parameterValues
        auto it = channelParameters.find(channel);
        if (it != channelParameters.end()) {
            for (const auto& param : it->second) {
	      parameterValues[param] = 1.; //initialize
	      cout << "parameter " << param << " added with value: " << parameterValues[param] << endl;
            }
        } else {
            std::cerr << "Channel " << channel << " not found in channelParameters." << std::endl;
        }
    }

    cout << "Declare Parameters called correctly" << endl;
    return parameterValues;
}

Parameter goldenmodes::getPar(const std::string& baseName) const {
    // Look for the real and imaginary parts in the parameterValues map
    auto it_real = parameterValues.find(baseName + "_re");
    auto it_imag = parameterValues.find(baseName + "_im");

    if (it_real != parameterValues.end() && it_imag != parameterValues.end()) {
        // Construct a complex number from the real and imaginary parts
        return std::complex<double>(it_real->second, it_imag->second);
    } else {
        throw std::runtime_error("Error: Real or imaginary part for parameter " + baseName + " not found.");
    }
}

// Setter function: sets the value for a given parameter in the map
void goldenmodes::SetParameterValue(const std::string& paramName, double value) {
    parameterValues[paramName] = value;  // Insert or update the value for the given parameter name
}


//---------------------------------------------------------------------

double goldenmodes::getParameterValue(const std::string& paramName) const {
    // Find the parameter in the parameterValues map
    auto it = parameterValues.find(paramName);

    // If the parameter is found, return its value
    if (it != parameterValues.end()) {
        return it->second;
    } else {
        // If the parameter is not found, throw an error or return a default value
        std::cerr << "Error: Parameter " << paramName << " not found in parameterValues map." << std::endl;
        throw std::runtime_error("Parameter not found: " + paramName);
    }
}

//----------------------------------------------------------

 //compute decay amplitude for each channel
 void goldenmodes::compute_decay_amplitudes(const std::string& channel, bool conjugate) {
   //cout << "computing decay amplitude for channel " << channel << endl;
   amplitudes[channel] = std::complex<double>(0.0, 0.0);

    Parameter amp;
    Parameter amp_0;
    Parameter amp_paral;
    Parameter amp_perp;
   // Get the CKM elements for the current channel, apply conjugation based on the 'conjugate' flag
    std::complex<double> Vcd = conjugate ? ckm.getVdc() : ckm.getVcd();
    std::complex<double> Vcs = conjugate ? ckm.getVsc() : ckm.getVcs();
    std::complex<double> Vbc = conjugate ? ckm.getVcb() : ckm.getVbc();
    std::complex<double> Vtd = conjugate ? ckm.getVdt() : ckm.getVtd();
    std::complex<double> Vbt = conjugate ? ckm.getVtb() : ckm.getVbt();
    std::complex<double> Vud = conjugate ? ckm.getVdu() : ckm.getVud();
    std::complex<double> Vbu = conjugate ? ckm.getVub() : ckm.getVbu();
    std::complex<double> Vts = conjugate ? ckm.getVst() : ckm.getVts();
    std::complex<double> Vus = conjugate ? ckm.getVsu() : ckm.getVus();


    if (channel == "Bdjpsik0") {
      amp = Vcd * Vbc * getPar("A_Bdjpsik0")
	+ Vud*Vbu * getPar("B_Bdjpsik0");
        amplitudes[channel] = amp;
    } else if (channel == "Bdjpsip0") {
      amp = Vcs*Vbc*getPar("A_Bdjpsip0")
	+ Vus*Vbu*getPar("B_Bdjpsip0");
        amplitudes[channel] = amp;
    } else if (channel == "Bdjpsiom") {
      amp_0 = Vcd*Vbc*getPar("A_Bdjpsiom_0")
	+ Vud*Vbu*getPar("B_Bdjpsiom_0");
  amp_paral = Vcd*Vbc*getPar("A_Bdjpsiom_paral")
  + Vud*Vbu*getPar("B_Bdjpsiom_paral");
  amp_perp = Vcd*Vbc*getPar("A_Bdjpsiom_perp")
+ Vud*Vbu*getPar("B_Bdjpsiom_perp");
      amp = amp_0 + amp_perp + amp_paral;
        amplitudes[channel] = amp;
        amplitudes[channel + "_0"] = amp_0;
        amplitudes[channel + "_paral"] = amp_paral;
        amplitudes[channel + "_perp"] = amp_perp;
    } else if (channel == "Bpjpsikp") {
	  amp = Vcs * Vbc * getPar("A_Bpjpsikp")
	    + Vus*Vbu * getPar("B_Bpjpsikp");
        amplitudes[channel] = amp;
    } else if (channel == "Bpjpsipp") {
	  amp = Vcd*Vbc*getPar("A_Bpjpsipp")
	+ Vud*Vbu*getPar("B_Bpjpsipp");
        amplitudes[channel] = amp;
    } else if (channel == "Bsjpsiph") {
      amp_0 = Vcs*Vbc*getPar("A_Bsjpsiph_0")
+ Vus*Vbu*getPar("B_Bsjpsiph_0");
amp_paral = Vcs*Vbc*getPar("A_Bsjpsiph_paral")
+ Vus*Vbu*getPar("B_Bsjpsiph_paral");
amp_perp = Vcs*Vbc*getPar("A_Bsjpsiph_perp")
+ Vus*Vbu*getPar("B_Bsjpsiph_perp");
amp = amp_0 + amp_perp + amp_paral;
amplitudes[channel] = amp;
amplitudes[channel + "_0"] = amp_0;
amplitudes[channel + "_paral"] = amp_paral;
amplitudes[channel + "_perp"] = amp_perp;
    } else if (channel == "Bsjpsik0") {
	  amp = Vcd*Vbc*getPar("A_Bsjpsik0")
	+ Vud*Vbu*getPar("B_Bsjpsik0");
        amplitudes[channel] = amp;
    } else if (channel == "Bdjpsikst") {
      amp_0 = Vcd*Vbc*getPar("A_Bdjpsikst_0")
	+ Vud*Vbu*getPar("B_Bdjpsikst_0");
  amp_paral = Vcd*Vbc*getPar("A_Bdjpsikst_paral")
  + Vud*Vbu*getPar("B_Bdjpsikst_paral");
  amp_perp = Vcd*Vbc*getPar("A_Bdjpsikst_perp")
+ Vud*Vbu*getPar("B_Bdjpsikst_perp");
      amp = amp_0 + amp_perp + amp_paral;
      amplitudes[channel] = amp;
      amplitudes[channel + "_0"] = amp_0;
      amplitudes[channel + "_paral"] = amp_paral;
      amplitudes[channel + "_perp"] = amp_perp;
    } else if (channel == "Bdjpsirh") {
      amp_0 = Vcd*Vbc*getPar("A_Bdjpsirh_0")
	+ Vud*Vbu*getPar("B_Bdjpsirh_0");
  amp_paral = Vcd*Vbc*getPar("A_Bdjpsirh_paral")
  + Vud*Vbu*getPar("B_Bdjpsirh_paral");
  amp_perp = Vcd*Vbc*getPar("A_Bdjpsirh_perp")
+ Vud*Vbu*getPar("B_Bdjpsirh_perp");
      amp = amp_0 + amp_perp + amp_paral;
      amplitudes[channel] = amp;
      amplitudes[channel + "_0"] = amp_0;
      amplitudes[channel + "_paral"] = amp_paral;
      amplitudes[channel + "_perp"] = amp_perp;
    } else if (channel == "Bsjpsikst") {
      amp_0 = Vcs*Vbc*getPar("A_Bsjpsikst_0")
+ Vus*Vbu*getPar("B_Bsjpsikst_0");
amp_paral = Vcs*Vbc*getPar("A_Bsjpsikst_paral")
+ Vus*Vbu*getPar("B_Bsjpsikst_paral");
amp_perp = Vcs*Vbc*getPar("A_Bsjpsikst_perp")
+ Vus*Vbu*getPar("B_Bsjpsikst_perp");
amp = amp_0 + amp_perp + amp_paral;
amplitudes[channel] = amp;
amplitudes[channel + "_0"] = amp_0;
amplitudes[channel + "_paral"] = amp_paral;
amplitudes[channel + "_perp"] = amp_perp;
} else {
      cout << "WARNING: amplitude for channel " << channel << " not found" << endl;
    }
 }

 // Getter for amplitudes
 Parameter goldenmodes::get_amplitude(const std::string& channel) {
     // First, check if the amplitude is already computed
     auto it = amplitudes.find(channel);
     if (it != amplitudes.end()) {
         return it->second;
     }

     // Extract the base channel name (if it's a polarized component)
     std::string base_channel = channel;
     if (channel.find("_0") != std::string::npos) {
         base_channel = channel.substr(0, channel.find("_0"));
     } else if (channel.find("_paral") != std::string::npos) {
         base_channel = channel.substr(0, channel.find("_paral"));
     } else if (channel.find("_perp") != std::string::npos) {
         base_channel = channel.substr(0, channel.find("_perp"));
     }

     // Compute the decay amplitude for the base channel if not already computed
     if (amplitudes.find(base_channel) == amplitudes.end()) {
         compute_decay_amplitudes(base_channel, false);
     }

     // Check if the specific requested amplitude exists
     it = amplitudes.find(channel);
     if (it == amplitudes.end()) {
         std::cerr << "ERROR: Amplitude not found for channel: " << channel << std::endl;
         std::cerr << "Available amplitudes: ";
         for (const auto& pair : amplitudes) {
             std::cerr << pair.first << " ";
         }
         std::cerr << std::endl;
         throw std::runtime_error("Amplitude not found for channel: " + channel);
     }

     return it->second;
 }





 // Getter for conjugated amplitudes
Parameter goldenmodes::get_conjugate_amplitude(const std::string& channel) {
    // First, check if the conjugate amplitude already exists
    auto it = amplitudes.find(channel);
    if (it != amplitudes.end()) {
        return it->second;
    }

    // Extract the base channel name (if it's a polarized component)
    std::string base_channel = channel;
    if (channel.find("_0") != std::string::npos) {
        base_channel = channel.substr(0, channel.find("_0"));
    } else if (channel.find("_paral") != std::string::npos) {
        base_channel = channel.substr(0, channel.find("_paral"));
    } else if (channel.find("_perp") != std::string::npos) {
        base_channel = channel.substr(0, channel.find("_perp"));
    }

    // Compute the conjugate decay amplitude for the base channel if not already computed
    if (amplitudes.find(base_channel) == amplitudes.end()) {
        compute_decay_amplitudes(base_channel, true);
    }

    // Check if the specific conjugate amplitude exists
    it = amplitudes.find(channel);
    if (it == amplitudes.end()) {
        std::cerr << "ERROR: Conjugate amplitude not found for channel: " << channel << std::endl;
        std::cerr << "Available amplitudes: ";
        for (const auto& pair : amplitudes) {
            std::cerr << pair.first << " ";
        }
        std::cerr << std::endl;
        throw std::runtime_error("Conjugate amplitude not found for channel: " + channel);
    }

    return it->second;
}



// Destructor
goldenmodes::~goldenmodes() {
    // Ensure any dynamically allocated resources are properly cleaned
    // delete histos;  // Uncomment if histos is dynamically allocated
}

//-----------------------------------------------------------
// Helper function to parse a decay channel name
std::pair<std::string, std::pair<std::string, std::string>>
goldenmodes::parseChannel(const std::string& channel) const {
    // Identify the B meson part
    std::string bMeson;
    if (channel.rfind("Bp", 0) == 0) {
        bMeson = "Bp";
    } else if (channel.rfind("Bd", 0) == 0) {
        bMeson = "Bd";
    } else if (channel.rfind("Bs", 0) == 0) {
        bMeson = "Bs";
    } else {
        throw std::runtime_error("Error in parseChannel: Unknown B meson in channel: " + channel);
    }

    // The remaining part should start with "jpsi" (or a variation)
    const std::string jpsiIdentifier = "jpsi";
    std::string remaining = channel.substr(bMeson.length());

    if (remaining.rfind(jpsiIdentifier, 0) != 0) {
        throw std::runtime_error("Error in parseChannel: Expected 'jpsi' in channel: " + channel);
    }

    // Extract the second final-state meson
    remaining = remaining.substr(jpsiIdentifier.length());
    std::string meson2;

    // Look for a valid meson name in mesonMasses
    for (const auto& meson : mesonMasses) {
        if (remaining.rfind(meson.first, 0) == 0) {
            meson2 = meson.first;
            break;
        }
    }

    if (meson2.empty()) {
        throw std::runtime_error("Error in parseChannel: Unable to parse second final-state meson in channel: " + channel);
    }

    return {bMeson, {jpsiIdentifier, meson2}};
}

//----------------------------------------------------------
// Getter for B meson lifetime
double goldenmodes::getBMesonLifetime(const std::string& bMeson) const {
    static const std::unordered_map<std::string, double> lifetimes = {
        {"Bp", tau_Bp},
        {"Bd", tau_Bd},
        {"Bs", tau_Bs}
    };

    auto it = lifetimes.find(bMeson);
    if (it != lifetimes.end()) {
        return it->second;
    }

    throw std::runtime_error("Error in getBMesonLifetime: Unknown B meson '" + bMeson + "'");
}

//----------------------------------------------------------
// Getter for B meson mass
double goldenmodes::getBMesonMass(const std::string& bMeson) const {
    static const std::unordered_map<std::string, double> masses = {
        {"Bp", m_Bp},
        {"Bd", m_Bd},
        {"Bs", m_Bs}
    };

    auto it = masses.find(bMeson);
    if (it != masses.end()) {
        return it->second;
    }

    throw std::runtime_error("Error in getBMesonMass: Unknown B meson '" + bMeson + "'");
}


double goldenmodes::CalculateBR(Parameter amplitude, const std::string& channel) const {
    // Parse the channel to extract meson components
    auto parsed = parseChannel(channel);
    const std::string& bMeson = parsed.first;
    const std::string& meson1 = parsed.second.first;
    const std::string& meson2 = parsed.second.second;

    // Get the masses of the decaying B meson and final-state mesons
    double m_B  = getBMesonMass(bMeson);
    double m1   = mesonMasses.at(meson1);
    double m2   = mesonMasses.at(meson2);

    // Ensure physical kinematics: prevent sqrt of a negative number
    double mass_term1 = (m_B * m_B - (m1 + m2) * (m1 + m2));
    double mass_term2 = (m_B * m_B - (m1 - m2) * (m1 - m2));

    if (mass_term1 < 0 || mass_term2 < 0) {
        throw std::runtime_error("Error in CalculateBR: Kinematically forbidden decay for channel " + channel);
    }

    // Compute the magnitude of the final-state momentum
    double p = sqrt(mass_term1 * mass_term2) / (2.0 * m_B);

    // Compute the decay width (Γ)
    double decay_width = (std::pow(G_F, 2) / (32.0 * M_PI * h_t)) * (p / std::pow(m_B, 2)) * std::norm(amplitude);

    // Compute the branching ratio using the B meson's lifetime
    double lifetime = getBMesonLifetime(bMeson);
    double BR = decay_width * lifetime;

    // Apply symmetry factor for identical particles in the final state
    if (meson1 == meson2) {
        BR *= 0.5;
    }

    return BR;
}

// -----------------------------------------------------------------------
// Function to calculate A_CP asymmetry
double goldenmodes::CalculateAcp(const Parameter& amplitude, const Parameter& conjugate_amplitude) const {
    // Ensure amplitudes are non-zero to avoid division errors
    double A2 = std::norm(amplitude);
    double Abar2 = std::norm(conjugate_amplitude);

    if (A2 + Abar2 == 0) {
        throw std::runtime_error("Error: CalculateAcp - Both amplitudes are zero, division by zero detected.");
    }

    // Compute A_CP
    return (Abar2 - A2) / (A2 + Abar2);
}

// ---------------------------------------------------------------------
// Function to calculate direct CP violation parameter C
double goldenmodes::CalculateC(const Parameter& amplitude, const Parameter& conjugate_amplitude, const std::string& channel) {
    // Parse the channel to determine the B meson type
    auto parsed = parseChannel(channel);
    std::string bMeson = parsed.first;

    // Get q/p ratio for Bd or Bs
    std::complex<double> q_p = (bMeson == "Bd") ? ckm.get_q_p_Bd() : ckm.get_q_p_Bs();

    // Special case for K0s and K0l channels (apply q/p_KS)
    if (channel == "Bdjpsik0s" || channel == "Bdjpsik0l" || channel == "Bsjpsik0s") {
        q_p *= ckm.get_q_p_KS();  // Multiply by q/p for K0 mixing
    }

    // Get CP eigenvalue for the channel (ensure it exists)
    if (cpEigenvalue.find(channel) == cpEigenvalue.end()) {
        throw std::runtime_error("Error: CalculateC - CP eigenvalue missing for channel: " + channel);
    }
    double eta = cpEigenvalue.at(channel);

    // Compute λ (lambda) = η * (q/p) * (A_conjugate / A)
    if (std::abs(amplitude) == 0) {
        std::cerr << "Error: CalculateC - Zero amplitude for " << channel << ", division by zero detected." << std::endl;
        return 0.0; // Return safe value instead of crashing
    }

    std::complex<double> lambda = eta * q_p * (conjugate_amplitude / amplitude);

    // Compute C observable: C = (1 - |λ|^2) / (1 + |λ|^2)
    double mod_lambda_squared = std::norm(lambda);
    return (1.0 - mod_lambda_squared) / (1.0 + mod_lambda_squared);
}


// ---------------------------------------------------------------------
// Function to calculate CP violation parameter S
double goldenmodes::CalculateS(const Parameter& amplitude, const Parameter& conjugate_amplitude, const std::string& channel) {
    // Parse the channel to determine the B meson type
    auto parsed = parseChannel(channel);
    std::string bMeson = parsed.first;

    // Get q/p ratio for Bd or Bs
    std::complex<double> q_p = (bMeson == "Bd") ? ckm.get_q_p_Bd() : ckm.get_q_p_Bs();

    // Special case for K0s and K0l channels (apply q/p_KS)
    if (channel == "Bdjpsik0s" || channel == "Bdjpsik0l" || channel == "Bsjpsik0s") {
        q_p *= ckm.get_q_p_KS();  // Multiply by q/p for K0 mixing
    }

    // Get CP eigenvalue for the channel (ensure it exists)
    if (cpEigenvalue.find(channel) == cpEigenvalue.end()) {
        throw std::runtime_error("Error: CalculateS - CP eigenvalue missing for channel: " + channel);
    }
    double eta = cpEigenvalue.at(channel);

    // Compute λ (lambda) = η * (q/p) * (A_conjugate / A)
    if (std::abs(amplitude) == 0) {
        std::cerr << "Error: CalculateS - Zero amplitude for " << channel << ", division by zero detected." << std::endl;
        return 0.0; // Return safe value instead of crashing
    }

    std::complex<double> lambda = eta * q_p * (conjugate_amplitude / amplitude);

    // Compute S observable: S = 2 Im(λ) / (1 + |λ|^2)
    double mod_lambda_squared = std::norm(lambda);
    return (2.0 * std::imag(lambda)) / (1.0 + mod_lambda_squared);
}



    std::pair<double, double> goldenmodes::CalculatePhiAndLambda(const Parameter& amplitude, const Parameter& conjugate_amplitude, const std::string& channel) {
       // Ensure the amplitude is nonzero to avoid division by zero
       if (std::abs(amplitude) == 0) {
           throw std::runtime_error("Error: CalculatePhiAndLambda - Zero amplitude detected for channel: " + channel);
       }

       // Get q/p for the Bs meson
       std::complex<double> q_p = ckm.get_q_p_Bs();

       // Compute lambda = (q/p) * (A_cp / A_conj)
       std::complex<double> lambda = q_p * (conjugate_amplitude / amplitude);

       // Compute |lambda|
       double mod_lambda = std::abs(lambda);

       // Compute phi_s = -arg(lambda)
       double phi_s = -std::arg(lambda);

       return {phi_s, mod_lambda};
   }



std::pair<std::vector<std::string>, std::string> goldenmodes::extractChannelFromCorrKey(const std::string& corr_key) {
    std::vector<std::string> channels;
    std::string experiment;

    if (corr_key.rfind("CS_", 0) == 0) {
        // Format: "CS_channel_exp"
        size_t underscore = corr_key.find("_", 3);
        if (underscore != std::string::npos) {
            channels.push_back(corr_key.substr(3, underscore - 3));  // Extract channel
            experiment = corr_key.substr(underscore + 1);  // Extract experiment
        }
    }
    else if (corr_key.rfind("ACP_", 0) == 0) {
        // Format: "ACP_channel1_channel2_exp"
        size_t first_underscore = corr_key.find("_", 4);
        size_t second_underscore = corr_key.find("_", first_underscore + 1);
        if (first_underscore != std::string::npos && second_underscore != std::string::npos) {
            channels.push_back(corr_key.substr(4, first_underscore - 4));
            channels.push_back(corr_key.substr(first_underscore + 1, second_underscore - first_underscore - 1));
            experiment = corr_key.substr(second_underscore + 1);
        }
    }
    else if (corr_key.rfind("phi_lambda_", 0) == 0) {
        // Format: "phi_lambda_channel_exp"
        size_t underscore = corr_key.find("_", 11);
        if (underscore != std::string::npos) {
            channels.push_back(corr_key.substr(11, underscore - 11));
            experiment = corr_key.substr(underscore + 1);
        }
    }
    else if (corr_key.rfind("phi_", 0) == 0) {
        // Format: "phi_channel_exp"
        size_t underscore = corr_key.find("_", 4);
        if (underscore != std::string::npos) {
            channels.push_back(corr_key.substr(4, underscore - 4));
            experiment = corr_key.substr(underscore + 1);
        }
    }
    else if (corr_key.rfind("polarization_", 0) == 0) {
        // Format: "polarization_channel_exp"
        size_t underscore = corr_key.find("_", 13);
        if (underscore != std::string::npos) {
            channels.push_back(corr_key.substr(13, underscore - 13));
            experiment = corr_key.substr(underscore + 1);
        }
    }
    else {
        std::cerr << "Warning!! Unknown key format: " << corr_key << std::endl;
    }

    return {channels, experiment};
}

std::map<std::string, double> goldenmodes::getPolarizationParams(
    const std::string& channel,
    const std::map<std::string, std::pair<Parameter, Parameter>>& amplitude_map)
{
    std::map<std::string, double> polarization_pars;

    try {
        // Get parameters for phases
        std::complex<double> B0    = getPar("B_" + channel + "_0");
        std::complex<double> Bperp = getPar("B_" + channel + "_perp");
        std::complex<double> Bparal= getPar("B_" + channel + "_paral");

        double delta_0      = std::arg(B0);
        double deltaparal   = std::arg(Bparal);
        double deltaperp    = std::arg(Bperp);
        double delta_paral  = deltaparal - delta_0;
        double delta_perp   = deltaperp - delta_0;

        // Get amplitude information
        auto amp_it = amplitude_map.find(channel);
        if (amp_it == amplitude_map.end()) {
            throw std::runtime_error("Amplitude map does not contain channel: " + channel);
        }

        std::complex<double> amp       = amp_it->second.first;
        std::complex<double> conj_amp  = amp_it->second.second;
        std::complex<double> avg_amp   = (amp + conj_amp) * 0.5;
        double norm_amp           = std::norm(avg_amp);

        if (norm_amp == 0) {
            throw std::runtime_error("Normalization factor is zero for channel: " + channel);
        }

        // Get decay amplitudes
        std::complex<double> amp_0      = get_amplitude(channel + "_0");
        std::complex<double> amp_paral  = get_amplitude(channel + "_paral");
        std::complex<double> amp_perp   = get_amplitude(channel + "_perp");
        std::complex<double> conj_amp_0     = get_conjugate_amplitude(channel + "_0");
        std::complex<double> conj_amp_paral = get_conjugate_amplitude(channel + "_paral");
        std::complex<double> conj_amp_perp  = get_conjugate_amplitude(channel + "_perp");

        std::complex<double> avg_A0      = (amp_0 + conj_amp_0) * 0.5;
        std::complex<double> avg_Aparal  = (amp_paral + conj_amp_paral) * 0.5;
        std::complex<double> avg_Aperp   = (amp_perp + conj_amp_perp) * 0.5;

        double norm_A0      = std::norm(avg_A0);
        double norm_Aparal  = std::norm(avg_Aparal);
        double norm_Aperp   = std::norm(avg_Aperp);

        // Calculate polarization fractions
        double f_0      = norm_A0 / norm_amp;
        double f_perp   = norm_Aperp / norm_amp;
        double f_paral  = norm_Aparal / norm_amp;

        // Store the results in the map
        polarization_pars["f_0_" + channel]      = f_0;
        polarization_pars["f_perp_" + channel]   = f_perp;
        polarization_pars["f_paral_" + channel]  = f_paral;
        polarization_pars["delta_paral_" + channel] = delta_paral;
        polarization_pars["delta_perp_" + channel]  = delta_perp;

    } catch (const std::exception& e) {
        std::cerr << "Error in getPolarizationParams for channel " << channel << ": " << e.what() << std::endl;
    }

    return polarization_pars;
}


double goldenmodes::Calculate_UncorrelatedObservables(const std::map<std::string, std::pair<Parameter, Parameter>>& amplitude_map) {
    double ll_uncorr = 0.0;  // Initialize log-likelihood contribution

    std::vector<std::string> channels = {
        "Bpjpsikp", "Bpjpsipp", "Bdjpsip0", "Bdjpsik0s", "Bdjpsik0l", "Bdjpsikst", "Bsjpsiph"
    };

    for (const auto& channel : channels) {
        try {
            auto it = amplitude_map.find(channel);
            if (it == amplitude_map.end()) {
                std::cerr << "Warning: Amplitude not found for " << channel << std::endl;
                continue;
            }
            const auto& amp_pair = it->second;

            // **Compute ACP, C, and S only if they exist in meas**
            if (meas.find("ACP" + channel) != meas.end()) {
                double acp = CalculateAcp(amp_pair.first, amp_pair.second);
                obs["ACP_" + channel] = acp;

                double observed = meas.at("ACP" + channel).getMean();
                double uncertainty = meas.at("ACP" + channel).getSigma();
                double diff = acp - observed;
                ll_uncorr += -0.5 * (diff * diff / (uncertainty * uncertainty));
            }

            if (meas.find("C" + channel) != meas.end()) {
                double c = CalculateC(amp_pair.first, amp_pair.second, channel);
                obs["C_" + channel] = c;

                double observed = meas.at("C" + channel).getMean();
                double uncertainty = meas.at("C" + channel).getSigma();
                double diff = c - observed;
                ll_uncorr += -0.5 * (diff * diff / (uncertainty * uncertainty));
            }

            if (meas.find("S" + channel) != meas.end()) {
                double s = CalculateS(amp_pair.first, amp_pair.second, channel);
                obs["S_" + channel] = s;

                double observed = meas.at("S" + channel).getMean();
                double uncertainty = meas.at("S" + channel).getSigma();
                double diff = s - observed;
                ll_uncorr += -0.5 * (diff * diff / (uncertainty * uncertainty));
            }

            // **Compute and store branching ratios if they exist in meas**
            std::string br_obsKey = "BR_" + channel;
            std::string br_measKey = "BR" + channel;
            if (meas.find(br_measKey) != meas.end()) {
                double br_predicted = CalculateBR(amp_pair.first, channel);
                obs[br_obsKey] = br_predicted;

                double observed = meas.at(br_measKey).getMean();
                double uncertainty = meas.at(br_measKey).getSigma();
                double diff = br_predicted - observed;
                ll_uncorr += -0.5 * (diff * diff / (uncertainty * uncertainty));
            }

        } catch (const std::exception& e) {
            std::cerr << "Error in Calculate_UncorrelatedObservables for " << channel << ": " << e.what() << std::endl;
        }
    }

    // **Handle BR Ratios (`R_`)**
    std::vector<std::pair<std::string, std::pair<std::string, std::string>>> br_ratios = {
        {"R_Bpjpsipp_Bpjpsikp", {"Bpjpsipp", "Bpjpsikp"}},
        {"R_Bdjpsiom_Bdjpsirh", {"Bdjpsiom", "Bdjpsirh"}},
        {"R_Bdjpsikst_Bdjpsik0", {"Bdjpsikst", "Bdjpsik0"}}
    };

    for (const auto& [ratioKey, channels] : br_ratios) {
        if (obs.find("BR_" + channels.first) != obs.end() && obs.find("BR_" + channels.second) != obs.end()) {
            double BR1 = obs["BR_" + channels.first];
            double BR2 = obs["BR_" + channels.second];

            double R_predicted = BR1 / BR2;
            obs[ratioKey] = R_predicted;

            if (meas.find(ratioKey) != meas.end()) {
                double observed = meas.at(ratioKey).getMean();
                double uncertainty = meas.at(ratioKey).getSigma();
                double diff = R_predicted - observed;
                ll_uncorr += -0.5 * (diff * diff / (uncertainty * uncertainty));
            }
        } else {
            std::cerr << "Warning: BR missing for " << ratioKey << std::endl;
        }
    }

    // **Handle Delta A (`deltaA_`)**
    std::vector<std::pair<std::string, std::pair<std::string, std::string>>> deltaA_ratios = {
        {"deltaA_Bpjpsipp_Bpjpsikp", {"ACP_Bpjpsipp", "ACP_Bpjpsikp"}}
    };

    for (const auto& [deltaAKey, channels] : deltaA_ratios) {
        if (obs.find(channels.first) != obs.end() && obs.find(channels.second) != obs.end()) {
            double ACP1 = obs[channels.first];
            double ACP2 = obs[channels.second];

            double deltaA_predicted = ACP1 - ACP2;
            obs[deltaAKey] = deltaA_predicted;

            if (meas.find(deltaAKey) != meas.end()) {
                double observed = meas.at(deltaAKey).getMean();
                double uncertainty = meas.at(deltaAKey).getSigma();
                double diff = deltaA_predicted - observed;
                ll_uncorr += -0.5 * (diff * diff / (uncertainty * uncertainty));
            }
        } else {
            std::cerr << "Warning: ACP missing for " << deltaAKey << std::endl;
        }
    }

    return ll_uncorr;
}



 //----------------------------------------------------------------------------------

double goldenmodes::Calculate_CorrelatedObservables(const std::map<std::string, std::pair<Parameter, Parameter>>& amplitude_map) {
    double ll_corr = 0.0;  // Initialize the log-likelihood contribution
    TVectorD corr(2);  // For correlated observables (e.g., C and S)
    TVectorD corr4(4);
    TVectorD corr5(5);
    TVectorD corr6(6);


    // **Correlated observables for Bdjpsik0s**
    {
        const auto& amp_pair_Bdjpsik0s = amplitude_map.at("Bdjpsik0s");
        if (obs.find("C_Bdjpsik0s") != obs.end() && obs.find("S_Bdjpsik0s") != obs.end()) {
          corr(0) = obs["C_Bdjpsik0s"];
          corr(1) = obs["S_Bdjpsik0s"];

          ll_corr += corrmeas.at("CS_Bdjpsik0s_LHCb2023").logweight(corr);
          ll_corr += corrmeas.at("CS_Bdjpsik0s_BelleII2024").logweight(corr);
        } else {
          std::cerr << "Error: C and S values for Bdjpsik0s not found in obs map!" << std::endl;
        }

    }

    // **Correlated observables for Bdjpsip0**
    {
        const auto& amp_pair_Bdjpsip0 = amplitude_map.at("Bdjpsip0");
        if (obs.find("C_Bdjpsip0") != obs.end() && obs.find("S_Bdjpsip0") != obs.end()) {
          corr(0) = obs["C_Bdjpsip0"];
          corr(1) = obs["S_Bdjpsip0"];
          ll_corr += corrmeas.at("CS_Bdjpsip0_BaBar2008").logweight(corr);
        } else {
          std::cerr << "Error: C and S values for Bdjpsip0 not found in obs map!" << std::endl;
        }
    }

    // **Correlated observables for Bdjpsirh**
    {
      try {
        const auto& amp_pair_Bdjpsirh = amplitude_map.at("Bdjpsirh");
        double c_predicted = CalculateC(amp_pair_Bdjpsirh.first, amp_pair_Bdjpsirh.second, "Bdjpsirh");
        double s_predicted = CalculateS(amp_pair_Bdjpsirh.first, amp_pair_Bdjpsirh.second, "Bdjpsirh");
        obs["C_Bdjpsirh"] = c_predicted;
        obs["S_Bdjpsirh"] = s_predicted;
        corr(0) = c_predicted;
        corr(1) = s_predicted;
        ll_corr += corrmeas.at("CS_Bdjpsirh_LHCb2014").logweight(corr);
      } catch (const std::out_of_range& e) { // Catch exception for missing key
        std::cerr << "Error: Bdjpsirh not found in amplitude_map! " << e.what() << std::endl;
      } catch (const std::exception& e) {
         std::cerr << "Unexpected error in correlated observables for Bdjpsirh: " << e.what() << std::endl;
        }
      }


    // **Correlated observables for Bsjpsiph**
    {
      auto pol_params = getPolarizationParams("Bsjpsiph", amplitude_map);

        // Extract the polarization fractions and phases:
        double f_0       = pol_params.at("f_0_Bsjpsiph");
        double f_paral   = pol_params.at("f_paral_Bsjpsiph");
        double f_perp    = pol_params.at("f_perp_Bsjpsiph");
        double delta_paral = pol_params.at("delta_paral_Bsjpsiph");
        double delta_perp  = pol_params.at("delta_perp_Bsjpsiph");

        obs["f_0_Bsjpsiph"]       = f_0;
        obs["f_paral_Bsjpsiph"]   = f_paral;
        obs["f_perp_Bsjpsiph"]    = f_perp;
        obs["delta_paral_Bsjpsiph"] = delta_paral;
        obs["delta_perp_Bsjpsiph"]  = delta_perp;

        // --- Retrieve additional parameters (e.g. weak phase and |lambda|) ---
        double phi_Bsjpsiph    = obs["phi_Bsjpsiph"];
        double lambda_Bsjpsiph = obs["lambda_Bsjpsiph"];


    // Loop over the correlated measurements
    for (auto& [key, corrObs] : corrmeas) {
        if (key.find("Bsjpsiph") != std::string::npos) {
            // Compute phi and lambda
            const auto& amp_pair_Bsjpsiph = amplitude_map.at("Bsjpsiph");
            //auto phi_lambda_result = CalculatePhiAndLambda(amp_pair_Bsjpsiph.first, amp_pair_Bsjpsiph.second, "Bsjpsiph");
            double phi_Bsjpsiph = obs["phi_Bsjpsiph"];
            double lambda_Bsjpsiph = obs["lambda_Bsjpsiph"];
            // Handle each correlated measurement key
            if (key == "phi_lambda_Bsjpsiph_LHCb2023") {
                TVectorD corr6(6);  // 6 correlated observables
                corr6(0) = phi_Bsjpsiph;      // Weak phase
                corr6(1) = lambda_Bsjpsiph;   // |lambda|
                corr6(2) = f_perp;            // Polarization fraction
                corr6(3) = f_0;               // Polarization fraction
                corr6(4) = delta_perp;        // Strong phase
                corr6(5) = delta_paral;       // Strong phase
                ll_corr += corrObs.logweight(corr6);
            } else if (key == "phi_Bsjpsiph_ATLAS2020B") {
                TVectorD corr5(5);  // 5 correlated observables
                corr5(0) = phi_Bsjpsiph;
                corr5(1) = f_paral;
                corr5(2) = f_0;
                corr5(3) = delta_perp;
                corr5(4) = delta_paral;
                ll_corr += corrObs.logweight(corr5);
            } else if (key == "phi_lambda_Bsjpsiph_LHCb2021") {
                TVectorD corr6(6);  // 6 correlated observables
                corr6(0) = phi_Bsjpsiph;
                corr6(1) = lambda_Bsjpsiph;
                corr6(2) = f_perp;
                corr6(3) = f_0;
                corr6(4) = delta_perp;
                corr6(5) = delta_paral;
                ll_corr += corrObs.logweight(corr6);
            } else if (key == "phi_Bsjpsiph_CMS2020") {
                TVectorD corr5(5);  // 5 correlated observables
                corr5(0) = phi_Bsjpsiph;
                corr5(1) = f_perp;
                corr5(2) = f_0;
                corr5(3) = delta_perp;
                corr5(4) = delta_paral;
                ll_corr += corrObs.logweight(corr5);
            } else if (key == "phi_lambda_Bsjpsiph_LHCb2019" || key == "phi_lambda_Bsjpsiph_LHCb2017") {
                TVectorD corr2(2);  // Only phi and lambda
                corr2(0) = phi_Bsjpsiph;
                corr2(1) = lambda_Bsjpsiph;
                ll_corr += corrObs.logweight(corr2);
            } else if (key == "phi_Bsjpsiph_ATLAS2016") {
                TVectorD corr5(5);  // 5 correlated observables
                corr5(0) = phi_Bsjpsiph;
                corr5(1) = f_paral;
                corr5(2) = f_0;
                corr5(3) = delta_perp;
                corr5(4) = delta_paral;
                ll_corr += corrObs.logweight(corr5);
            } else if (key == "phi_Bsjpsiph_ATLAS2014") {
                TVectorD corr5(5);  // 5 correlated observables
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

  //correlated polarization measurments for C_Bdjpsikst
{
    auto pol_params = getPolarizationParams("Bdjpsikst", amplitude_map);

      // Extract the polarization fractions and phases:
      double f_0       = pol_params.at("f_0_Bdjpsikst");
      double f_paral   = pol_params.at("f_paral_Bdjpsikst");
      double f_perp    = pol_params.at("f_perp_Bdjpsikst");
      double delta_paral = pol_params.at("delta_paral_Bdjpsikst");
      double delta_perp  = pol_params.at("delta_perp_Bdjpsikst");

      obs["f_0_Bdjpsikst"]       = f_0;
      obs["f_paral_Bdjpsikst"]   = f_paral;
      obs["f_perp_Bdjpsikst"]    = f_perp;
      obs["delta_paral_Bdjpsikst"] = delta_paral;
      obs["delta_perp_Bdjpsikst"]  = delta_perp;

      TVectorD corr4(4);  // 4 correlated observables
      corr4(0) = f_paral;            // Polarization fraction
      corr4(1) = f_perp;               // Polarization fraction
      corr4(2) = delta_paral;        // Strong phase
      corr4(3) = delta_perp;       // Strong phase
      ll_corr += corrmeas.at("polarization_Bdjpsikst_LHCb2013").logweight(corr4);
    }


{
  auto pol_params = getPolarizationParams("Bsjpsikst", amplitude_map);

    // Extract the polarization fractions and phases:
    double f_0       = pol_params.at("f_0_Bsjpsikst");
    double f_paral   = pol_params.at("f_paral_Bsjpsikst");
    double delta_paral = pol_params.at("delta_paral_Bsjpsikst");
    double delta_perp  = pol_params.at("delta_perp_Bsjpsikst");

    obs["f_0_Bsjpsikst"]       = f_0;
    obs["f_paral_Bsjpsikst"]   = f_paral;
    obs["delta_paral_Bsjpsikst"] = delta_paral;
    obs["delta_perp_Bsjpsikst"]  = delta_perp;
  // --- Retrieve additional parameters  ---
  complex<double> amp_0      = get_amplitude("Bsjpsikst_0");
  complex<double> amp_paral  = get_amplitude("Bsjpsikst_paral");
  complex<double> amp_perp   = get_amplitude("Bsjpsikst_perp");
  complex<double> conj_amp_0     = get_conjugate_amplitude("Bsjpsikst_0");
  complex<double> conj_amp_paral = get_conjugate_amplitude("Bsjpsikst_paral");
  complex<double> conj_amp_perp  = get_conjugate_amplitude("Bsjpsikst_perp");


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


double goldenmodes::LogLikelihood(const std::vector<double>& parameters) {
    static int iteration_counter = 0;
    ++iteration_counter;

    obs.clear();  // Clear obs map for each iteration
    double ll = 0.0;  // Log-likelihood accumulator


    // Check that the parameter vector has the expected size
    int expectedSize = 0;
    for (const auto& channel : channelNamesSU3) {
        auto it = channelParameters.find(channel);
        if (it != channelParameters.end()) {
            expectedSize += it->second.size();
        }
    }
    // Last 4 parameters are for CKM
    if (parameters.size() != expectedSize + 4) {
        std::cerr << "Error: parameters.size() = " << parameters.size()
                  << ", but expected size = " << expectedSize + 4 << std::endl;
        return -1e30;
    }
    if (parameters.empty()) {
        std::cerr << "Error: Empty parameters vector!" << std::endl;
        return -1e30;
    }

    // Unpack CKM parameters and compute CKM elements
    std::vector<double> ckmParams(parameters.end() - 4, parameters.end());
    ckm.computeCKM(ckmParams[0], ckmParams[1], ckmParams[2], ckmParams[3], true);

    // Unpack the parameters for each channel (from channelNamesSU3)
    int index = 0;
    for (const auto& channel : channelNamesSU3) {
        auto it = channelParameters.find(channel);
        if (it != channelParameters.end()) {
            // Loop over parameter pairs (real and imaginary)
            for (size_t j = 0; j < it->second.size(); j += 2) {
                if (index + 1 >= expectedSize) { // Only loop over channel parameters
                    break;
                }

                double realPart = parameters[index];
                double imagPart = parameters[index + 1];
                const std::string& realParamName = it->second[j];
                const std::string& imagParamName = it->second[j + 1];

                SetParameterValue(realParamName, realPart);
                SetParameterValue(imagParamName, imagPart);

                index += 2;
            }
        } else {
            std::cerr << "Channel " << channel << " not found in channelParameters." << std::endl;
            return -1e30;
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
    std::map<std::string, std::pair<Parameter, Parameter>> amplitude_map;
    amplitude_map.clear();

    for (const std::string& channel : channelNamesSU3) {
        if (std::find(vectorMesonChannels.begin(), vectorMesonChannels.end(), channel) != vectorMesonChannels.end()) {
            // Vector meson channels: handle polarized parameters
            Parameter A_0 = getPar("A_" + channel + "_0");
            Parameter B_0 = getPar("B_" + channel + "_0");
            Parameter A_perp = getPar("A_" + channel + "_perp");
            Parameter B_perp = getPar("B_" + channel + "_perp");
            Parameter A_paral = getPar("A_" + channel + "_paral");
            Parameter B_paral = getPar("B_" + channel + "_paral");

            // Store real and imaginary parts
            obs["A_0_re_" + channel] = A_0.real();
            obs["A_0_im_" + channel] = A_0.imag();
            obs["B_0_re_" + channel] = B_0.real();
            obs["B_0_im_" + channel] = B_0.imag();

            obs["A_perp_re_" + channel] = A_perp.real();
            obs["A_perp_im_" + channel] = A_perp.imag();
            obs["B_perp_re_" + channel] = B_perp.real();
            obs["B_perp_im_" + channel] = B_perp.imag();

            obs["A_paral_re_" + channel] = A_paral.real();
            obs["A_paral_im_" + channel] = A_paral.imag();
            obs["B_paral_re_" + channel] = B_paral.real();
            obs["B_paral_im_" + channel] = B_paral.imag();
        } else {
            // Pseudoscalar channels: handle A and B
            Parameter A = getPar("A_" + channel);
            Parameter B = getPar("B_" + channel);
            obs["A_re_" + channel] = A.real();
            obs["A_im_" + channel] = A.imag();
            obs["B_re_" + channel] = B.real();
            obs["B_im_" + channel] = B.imag();
        }

        // Directly store original amplitudes
    amplitude_map[channel] = {get_amplitude(channel), get_conjugate_amplitude(channel)};

    // For CP observables, add K0s and K0l cases with the same amplitudes as Bdjpsik0
    if (channel == "Bdjpsik0") {
        amplitude_map["Bdjpsik0s"] = amplitude_map["Bdjpsik0"];
        amplitude_map["Bdjpsik0l"] = amplitude_map["Bdjpsik0"];
    } else if (channel == "Bsjpsik0") {
        amplitude_map["Bsjpsik0s"] = amplitude_map["Bsjpsik0"];
    }
}

    // Check NaN/Inf values in amplitudes
    for (const auto& [chan, amp_pair] : amplitude_map) {
        if (std::isnan(abs(amp_pair.first)) || std::isinf(abs(amp_pair.first)) ||
            std::isnan(abs(amp_pair.second)) || std::isinf(abs(amp_pair.second))) {
            std::cerr << "Invalid amplitude (NaN or Inf) for channel: " << chan << std::endl;
            return -100;
        }
    }

    // Compute Branching Ratios and Likelihood Contribution
    for (const std::string& channel : channelNames) {
        const auto& amp_pair = amplitude_map[channel];

        if (channel == "Bdjpsik0s" || channel == "Bdjpsik0l") {
            const auto& amp_pair = amplitude_map["Bdjpsik0"];  // Use the non-rotated state
            double br_A = CalculateBR(amp_pair.first, "Bdjpsik0");  // Compute BR with non-rotated states
            double br_conj_A = CalculateBR(amp_pair.second, "Bdjpsik0");
            double br_predicted = 0.5 * (br_A + br_conj_A);
            obs["BR_Bdjpsik0"] = br_predicted;

            try {
                std::string br_key = "BRBdjpsik0";
                double br_observed = meas.at(br_key).getMean();
                double br_uncertainty = meas.at(br_key).getSigma();
                double diff = br_predicted - br_observed;
                ll += -0.5 * (diff * diff / (br_uncertainty * br_uncertainty));
            } catch (const std::out_of_range&) {
                std::cerr << "Error: Branching ratio for Bdjpsik0 not found in meas map." << std::endl;
                return -1e30;
            }
        } else {
            // Calculate BR for each channel
            double br_A = CalculateBR(amp_pair.first, channel);
            double br_conj_A = CalculateBR(amp_pair.second, channel);
            double br_predicted = 0.5 * (br_A + br_conj_A);
            obs["BR_" + channel] = br_predicted;

            // Compute log-likelihood for BR
            try {
                std::string br_key = "BR" + channel;
                double br_observed = meas.at(br_key).getMean();
                double br_uncertainty = meas.at(br_key).getSigma();
                double diff = br_predicted - br_observed;
                ll += -0.5 * (diff * diff / (br_uncertainty * br_uncertainty));
            } catch (const std::out_of_range&) {
                std::cerr << "Error: Branching ratio for " << channel << " not found in meas map." << std::endl;
                return -1e30;
            }
        }
    }


    //Add contributions from uncorrelated and correlated observables
    ll += Calculate_UncorrelatedObservables(amplitude_map);
    ll += Calculate_CorrelatedObservables(amplitude_map);

    return ll;
}

//---------------------------------------------------------


void goldenmodes::MCMCUserIterationInterface() {
    // Loop over all MCMC chains
  const unsigned int log_interval = 1000;  // Log every 1000 iterations
  std::vector<double> pars;

    for (unsigned int i = 0; i < fMCMCNChains; ++i) {
        pars = fMCMCStates.at(i).parameters;
	try {
	  LogLikelihood(pars);
	} catch (const std::exception& e) {
	  std::cerr << "Error in LogLikelihood: " << e.what() << std::endl;
	  return;
	}

        // Directly fill histograms
        histos.fillh1d();
        histos.fillh2d();
    }
}



void goldenmodes::SaveHistograms(const std::string& filename) {
    TFile file(filename.c_str(), "RECREATE");

    // Write all histograms to the file
    for (const auto& histPair : histos.h1d) {
        histPair.second->Write();
    }

    // Write all 2D histograms
    for (const auto& histPair : histos.h2d) {
        histPair.second->Write();
    }

    file.Close();
    std::cout << "Histograms saved to " << filename << std::endl;
}

// ---------------------------------------------------------

void goldenmodes::PrintObservablePulls(const std::string& filename) {
    std::ofstream outfile(filename);

    if (!outfile.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return;
    }

    // Header for the output file
    outfile << "Observable\t\t Measurement\t\t Mean\t\t Sigma\t\t Pull\n";

    // Define mappings for observable names and measurement keys
    std::map<std::string, std::string> measMap = {
        {"BR_Bdjpsik0", "BRBdjpsik0"}, {"BR_Bdjpsip0", "BRBdjpsip0"}, {"BR_Bdjpsiom", "BRBdjpsiom"},
        {"BR_Bpjpsikp", "BRBpjpsikp"}, {"BR_Bpjpsipp", "BRBpjpsipp"}, {"BR_Bsjpsiph", "BRBsjpsiph"},
        {"BR_Bsjpsik0s", "BRBsjpsik0s"}, {"BR_Bdjpsirh", "BRBdjpsirh"}, {"BR_Bdjpsikst", "BRBdjpsikst"},
        {"BR_Bsjpsikst", "BRBsjpsikst"}, {"C_Bdjpsik0s", "CBdjpsik0s"}, {"C_Bdjpsik0l", "CBdjpsik0l"},
        {"S_Bdjpsik0s", "SBdjpsik0s"}, {"S_Bdjpsik0l", "SBdjpsik0l"}, {"C_Bdjpsip0", "CBdjpsip0"},
        {"S_Bdjpsip0", "SBdjpsip0"}, {"ACP_Bpjpsikp", "ACPBpjpsikp"}, {"ACP_Bpjpsipp", "ACPBpjpsipp"},
        {"deltaA_Bpjpsipp_Bpjpsikp", "deltaA_Bpjpsipp_Bpjpsikp"}, {"R_Bpjpsipp_Bpjpsikp", "R_Bpjpsipp_Bpjpsikp"},
        {"lambda_Bsjpsiph", "lambda_Bsjpsiph"}, {"phi_Bsjpsiph", "phi_Bsjpsiph"}, {"f_perp_Bsjpsiph", "f_perp_Bsjpsiph"},
        {"f_paral_Bsjpsiph", "f_paral_Bsjpsiph"}, {"f_0_Bsjpsiph", "f_0_Bsjpsiph"}, {"delta_perp_Bsjpsiph", "delta_perp_Bsjpsiph"},
        {"delta_paral_Bsjpsiph", "delta_paral_Bsjpsiph"}, {"R_Bdjpsiom_Bdjpsirh", "R_Bdjpsiom_Bdjpsirh"},
        {"R_Bdjpsikst_Bdjpsik0", "R_Bdjpsikst_Bdjpsik0"}, {"C_Bdjpsikst", "CBdjpsikst"}, {"S_Bdjpsikst", "SBdjpsikst"},
        {"f_0_Bdjpsikst", "f_0_Bdjpsikst"}, {"f_paral_Bdjpsikst", "f_paral_Bdjpsikst"}, {"f_perp_Bdjpsikst", "f_perp_Bdjpsikst"},
        {"delta_paral_Bdjpsikst", "delta_paral_Bdjpsikst"}, {"delta_perp_Bdjpsikst", "delta_perp_Bdjpsikst"},
        {"f_0_Bsjpsikst", "f_0_Bsjpsikst"}, {"f_paral_Bsjpsikst", "f_paral_Bsjpsikst"},
        {"delta_paral_Bsjpsikst", "delta_paral_Bsjpsikst"}, {"delta_perp_Bsjpsikst", "delta_perp_Bsjpsikst"},
        {"A0_CP_Bsjpsikst", "A0_CP_Bsjpsikst"}, {"Aparal_CP_Bsjpsikst", "Aparal_CP_Bsjpsikst"}, {"Aperp_CP_Bsjpsikst", "Aperp_CP_Bsjpsikst"}

    };

    for (const auto& histPair : histos.h1d) {
        const std::string& obs_name = histPair.first;
        TH1D* hist = histPair.second;

        double obs_mean = hist->GetMean();
        double obs_measurement = 0.0;
        double sigma_measurement = 0.0;
        bool found_measurement = false;

        // Check if the observable exists in `newmeas` first
        if (measMap.find(obs_name) != measMap.end()) {
            const std::string& key = measMap[obs_name];

            if (newmeas.find(key) != newmeas.end()) {
                obs_measurement = newmeas.at(key).getMean();
                sigma_measurement = newmeas.at(key).getSigma();
                found_measurement = true;
            } else if (meas.find(key) != meas.end()) {
                // Use `meas` only if it is not found in `newmeas`
                obs_measurement = meas.at(key).getMean();
                sigma_measurement = meas.at(key).getSigma();
                found_measurement = true;
            }
        }

        if (found_measurement) {
            // Calculate the pull value
            Pull pullCalculator(
                *hist,                   // Reference to the histogram
                hist->GetNbinsX(),       // Number of bins along X
                500,                     // Number of bins along Y
                0.0,                     // X-axis lower bound, triggers auto setting in makeCompatPlot
                0.0,                     // X-axis upper bound, triggers auto setting in makeCompatPlot
                0.0,                     // Y-axis lower bound, triggers auto setting in makeCompatPlot
                0.0,                     // Y-axis upper bound, triggers auto setting in makeCompatPlot
                true                     // Whether to use low stats or not
            );

            double pull = pullCalculator.calcPull(obs_measurement, sigma_measurement, true);
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
