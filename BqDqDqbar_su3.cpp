#include "BqDqDqbar.h"
#include "dato.h"
#include "PDGAverage.h"
#include "CKM.h"
#include "Pull.h"
#include "BaseMacros.h"

#include <TFile.h>

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

BqDqDqbar::BqDqDqbar() : BCModel() , histos(obs)
{
  cout << "constructor for BqDqDqbar_su(3) called: inserting the experimental data" << endl;
    std::vector<dato> data;
    PDGAverage pdgaverage;



    DeclareParameters();  // Ensure parameters are defined

      // Add CKM parameters directly in the constructor
    AddParameter("CKM_Vud", 0.97415, 0.97447);         // Vud parameter
    AddParameter("CKM_Vcb", 0.04046, 0.04194 );    // Vcb parameter
    AddParameter("CKM_Vub", 0.00349 , 0.00419);  // Vub parameter
    AddParameter("CKM_gamma", 1.12224671, 1.22347581);   // gamma parameter (in rad)

    SetPriorConstantAll();


    //measurments

    //BRBddpdm
    data.push_back(dato(2.12e-4, 0.16e-4, 0.18e-4)); //Belle:2012mef
    data.push_back(dato(2.8e-4, 0.4e-4, 0.5e-4)); //BaBar:2006uih

    pdgaverage.setData(data);
    pdgaverage.setName("BRBddpdm");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();
    if (meas.find("BRBddpdm") != meas.end()) {
    std::cout << "BRBddpdm inserted successfully" << std::endl;
}


    //BRBddpsdm
    data.push_back(dato(7.3e-3, 0.4e-3, 0.7e-3)); //Zupanc:2007pu
    data.push_back(dato(6.6e-3, 1.4e-3, 0.6e-3)); //BaBar:2006jvx
    data.push_back(dato(6.8e-3, 2.4e-3, 0.6e-3)); //CLEO:1995psi
    data.push_back(dato(10e-3, 9e-3, 1e-3)); //ARGUS:1991xej
    data.push_back(dato(5.3e-3, 3.0e-3, 0.5e-3)); //CLEO:1991roe


    pdgaverage.setData(data);
    pdgaverage.setName("BRBddpsdm");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();


    //BRBpdpd0b
    data.push_back(dato(3.85e-4, 0.31e-4, 0.38e-4)); //Belle:2008doh
    data.push_back(dato(3.8e-4, 0.6e-4, 0.5e-4)); //BaBar:2006uih

    pdgaverage.setData(data);
    pdgaverage.setName("BRBpdpd0b");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    //BRBpdpsd0b
    data.push_back(dato(8.6e-3, 0.2e-3, 1.1e-3));  //LHCb:2013sad
    data.push_back(dato(9.5e-3, 2.0e-3, 0.8e-3)); //BaBar:2006jvx
    data.push_back(dato(9.8e-3, 2.6e-3, 0.9e-3)); //CLEO:1995psi
    data.push_back(dato(14e-3, 8e-3, 1e-3)); //ARGUS:1991xej
    data.push_back(dato(13e-3, 6e-3, 1e-3)); //CLEO:1990mqz

    pdgaverage.setData(data);
    pdgaverage.setName("BRBpdpsd0b");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    //BRBsdpsdms
    data.push_back(dato(4.0e-3, 0.2e-3, 0.5e-3));  //LHCb:2013sad
    data.push_back(dato(5.9e-3, 1.0e-3, 1.3e-3)); //Belle:2012tsw
    data.push_back(dato(5.4e-3, 0.8e-3, 0.8e-3)); //CDF:2012xmd

    pdgaverage.setData(data);
    pdgaverage.setName("BRBsdpsdms");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();


    //BRBsdpdms
    meas.insert(pair<string, dato>("BRBsdpdms",dato(2.8e-4, 0.4e-4, 0.3e-4))); //LHCb:2014scu


    //ACP measurments

    //CBddpdm
    meas.insert(pair<string, dato>("CBddpdm",dato(-0.43, 0.16, 0.05))); //Belle:2012mef

    //SBddpdm
    meas.insert(pair<string, dato>("SBddpdm",dato(-0.99, 0.175, 0.08))); //Belle:2012mef

    //ACPBpdpd0b
    data.push_back(dato(0.00, 0.08, 0.02)); //Belle:2008doh
    data.push_back(dato(-0.13, 0.14, 0.02));  //BaBar:2006uih


    pdgaverage.setData(data);
    pdgaverage.setName("ACPBpdpd0b");
    pdgaverage.CalculateAverage();

    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();


    //mod(lambda) inserted in meas to predict phi_s
    meas.insert(pair<string, dato>("L_Bsdpsdms", dato(1.054, 0.099, 0.020))); //LHCb:2024gkk




    //------------------------------------------------
    std::vector<dato> CorrData;
    TMatrixDSym CorrStat(2);  // 2 measurements (C and S)
    TMatrixDSym CorrSyst(2);  // 2 measurements (C and S)


    //Bddpdm : C and S observables from BaBar:2008xnt
    CorrData.push_back(dato(-0.07, 0.23, 0.03));  // C observable
    CorrData.push_back(dato(-0.63, 0.36, 0.05)); // S observable


// Populate the correlation matrix
    CorrStat(0, 0) = 1.;       // Variance for C
    CorrStat(1, 1) = 1.;       // Variance for S
    CorrStat(0, 1) = -0.012;      // Correlation between C and S
    CorrStat(1, 0) = -0.012;      // Symmetric part (same as Corr(0, 1))


    //CorrSyst
    CorrSyst(0, 0) = 1.;       //
    CorrSyst(1, 1) = 1.;
    CorrSyst(0, 1) = 0.;      // Correlation between C and S
    CorrSyst(1, 0) = 0.;      // Symmetric part (same as Corr(0, 1))



// Insert correlated data for BaBar2008 into corrmeas
    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("CS_Bddpdm_BaBar2008",CorrelatedGaussianObservables(CorrData, CorrStat, CorrSyst)));

    CorrData.clear();
    //CorrStat.Clear();
    //CorrSyst.Clear();


    //Bddpdm : C and S observables from LHCb:2024gkk
    CorrData.push_back(dato(0.162, 0.088, 0.009));  // C observable
    CorrData.push_back(dato(-0.549, 0.085, 0.015)); // S observable


// Populate the correlation matrix
    CorrStat(0, 0) = 1.;       //
    CorrStat(1, 1) = 1.;
    CorrStat(0, 1) = 0.475;      // Correlation between C and S
    CorrStat(1, 0) = 0.475;      // Symmetric part (same as Corr(0, 1))


    // Populate the correlation matrix
    CorrSyst(0, 0) = 1.;       //
    CorrSyst(1, 1) = 1.;
    CorrSyst(0, 1) = 0.;      // Correlation between C and S
    CorrSyst(1, 0) = 0.;      // Symmetric part (same as Corr(0, 1))


// Insert correlated data for LHCb:2024gkk into corrmeas
    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("CS_Bddpdm_LHCb2024",CorrelatedGaussianObservables(CorrData, CorrStat, CorrSyst)));

    CorrData.clear();


    //ACPBpdpd0b and ACPBpdpsd0b from LHCb:2023wbb
    dato Bpdpd0b(2.5e-2, 1.0e-2, 0.4e-2, 0.3e-2);
    dato Bpdpsd0b(0.5e-2, 0.2e-2, 0.5e-2, 0.3e-2);

    CorrData.push_back(dato(Bpdpd0b.getMean(), Bpdpd0b.getSigma()));  // Bpdpd0b
    CorrData.push_back(dato(Bpdpsd0b.getMean(), Bpdpsd0b.getSigma())); // Bpdpsd0b

    TMatrixDSym Corr(2); //total correlation


// Populate the correlation matrix
    Corr(0, 0) = 1.0;       // Variance
    Corr(1, 1) = 1.0;       // Variance
    Corr(0, 1) = 0.386;      // Correlation
    Corr(1, 0) = 0.386;      // Symmetric part (same as Corr(0, 1))


// Insert correlated data for LHCb:2016inx into corrmeas
    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("ACP_Bpdpd0b_Bpdpsd0b_LHCb2023",CorrelatedGaussianObservables(CorrData, Corr)));

    CorrData.clear();

    data.push_back(dato(-0.43, 0.16, 0.05)); //Belle:2012mef
    data.push_back(dato(0.162, 0.088, 0.009)); //LHCb:2024gkk
    data.push_back(dato(-0.07, 0.23, 0.03));

    pdgaverage.setData(data);
    pdgaverage.setName("CBddpdm");
    pdgaverage.CalculateAverage();

    newmeas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    data.push_back(dato(-0.99, 0.175, 0.08));
    data.push_back(dato(-0.549, 0.085, 0.015));
    data.push_back(dato(-0.63, 0.36, 0.05));

    pdgaverage.setData(data);
    pdgaverage.setName("SBddpdm");
    pdgaverage.CalculateAverage();

    newmeas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    data.push_back(dato(0.00, 0.08, 0.02));
    data.push_back(dato(-0.13, 0.14, 0.02));
    data.push_back(dato(2.5e-2, 1.0e-2, 0.4e-2, 0.3e-2));

    pdgaverage.setData(data);
    pdgaverage.setName("ACPBpdpd0b");
    pdgaverage.CalculateAverage();

    newmeas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    newmeas.insert(pair<string, dato>("ACPBpdpsd0b", dato(0.5e-2, 0.2e-2, 0.5e-2, 0.3e-2)));

    std::cout << "All meas inserted" << endl;


   // Create histograms for uncorrelated observables (from meas)
    for (const auto& channel : channelNames) {
        // Branching ratio (BR)
        if (meas.find("BR" + channel) != meas.end()) {
            if (histos.h1d.find("BR_" + channel) == histos.h1d.end()) {  // Check if it already exists
                histos.createH1D("BR_" + channel, 500, 0.0, 0.0);

            }
        }

        // ACP, C, and S observables
        if (meas.find("ACP" + channel) != meas.end()) {
            if (histos.h1d.find("ACP_" + channel) == histos.h1d.end()) {
                histos.createH1D("ACP_" + channel, 500, 0.0, 0.0);

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

    histos.createH1D("phi_Bsdpsdms", 500, 0.0, 0.0);
    histos.createH1D("lamda_Bsdpsdms", 500, 0.0, 0.0);


    //create hinstogram to store delta_phi_s


    histos.createH1D("delta_phi", 500, 0.0, 0.0);
    // Create histograms for uncorrelated observables (from meas)

    for (const auto& channel : channelNames) {
      histos.createH1D("A_" + channel + "_re", 500, 0.0, 0.0);
      histos.createH1D("A_" + channel + "_im", 500, 0.0, 0.0);
      histos.createH1D("B_" + channel + "_re", 500, 0.0, 0.0);
      histos.createH1D("B_" + channel + "_im", 500, 0.0, 0.0);
      histos.createH1D("|A|(" + channel + ")", 500, 0.0, 0.0);
      histos.createH1D("A_phase (" + channel + ")", 500, 0.0, 0.0);
      histos.createH1D("|B|(" + channel + ")", 500, 0.0, 0.0);
      histos.createH1D("B_phase (" + channel +  ")", 500, 0.0, 0.0);
    }



    histos.createH1D("|E1|", 500, 0.0, 0.0);
    histos.createH1D("|A2|", 500, 0.0, 0.0);
    histos.createH1D("|A1|", 500, 0.0, 0.0);
    histos.createH1D("|P1_GIM|", 500, 0.0, 0.0);
    histos.createH1D("|P3_GIM|", 500, 0.0, 0.0);
    histos.createH1D("phase A2", 500, 0.0, 0.0);
    histos.createH1D("phase A1", 500, 0.0, 0.0);
    histos.createH1D("phase P1_GIM", 500, 0.0, 0.0);
    histos.createH1D("phase P3_GIM", 500, 0.0, 0.0);
    histos.createH1D("E1", 500, 0.0, 0.0);
    histos.createH1D("A2_re", 500, 0.0, 0.0);
    histos.createH1D("A2_im", 500, 0., 0.);
    histos.createH1D("A1_re", 500, 0.0, 0.0);
    histos.createH1D("A1_im", 500, 0., 0.);
    histos.createH1D("P1_GIM_re", 500, 0.0, 0.0);
    histos.createH1D("P1_GIM_im", 500, 0., 0.);
    histos.createH1D("P3_GIM_re", 500, 0.0, 0.0);
    histos.createH1D("P3_GIM_im", 500, 0., 0.);




}



//---------------------------------------------------------------


   //method to define the parameters needed to calculate each decay amplitude, uses BCModel method AddParameter
void BqDqDqbar::DefineParameters(const string& channel) {

  if (channel == "Bddpdm") {
        std::vector<std::string> params = {
            "E1_dcc_Bddpdm_re",
	    "E1_dcc_Bddpdm_im",
            "A2_cdc_Bddmdp_re",
	    "A2_cdc_Bddmdp_im",
            "P1_GIM_dc_Bddpdm_re",
	    "P1_GIM_dc_Bddpdm_im",
            "P3_GIM_dc_Bddmdp_re",
	    "P3_GIM_dc_Bddmdp_im"
        };
        channelParameters[channel] = params;
        for (const auto& param : params) {
            AddParameter(param, -40., 40.);
        }
	GetParameter("E1_dcc_Bddpdm_im").SetUpperLimit(0.0);
	GetParameter("E1_dcc_Bddpdm_im").SetLowerLimit(0.);
	GetParameter("E1_dcc_Bddpdm_re").SetLowerLimit(0.);
	GetParameter("P1_GIM_dc_Bddpdm_re").SetUpperLimit(130.);
	GetParameter("P1_GIM_dc_Bddpdm_re").SetLowerLimit(-130.);
	GetParameter("P3_GIM_dc_Bddmdp_re").SetUpperLimit(130.);
	GetParameter("P3_GIM_dc_Bddmdp_re").SetLowerLimit(-130.);
	GetParameter("P1_GIM_dc_Bddpdm_im").SetUpperLimit(100.);
	GetParameter("P1_GIM_dc_Bddpdm_im").SetLowerLimit(-100.);
	GetParameter("P3_GIM_dc_Bddmdp_im").SetUpperLimit(100.);
	GetParameter("P3_GIM_dc_Bddmdp_im").SetLowerLimit(-100.);


    }
  else if (channel == "Bddpsdm") {
        std::vector<std::string> params = {
            "su3_E1_scc_Bddpsdm_re",
	    "su3_E1_scc_Bddpsdm_im",
            "su3_P1_GIM_sc_Bddpsdm_re",
	    "su3_P1_GIM_sc_Bddpsdm_im"
	};
        channelParameters[channel] = params;
        for (const auto& param : params) {
	  AddParameter(param, -1., 1.);
        }
    } else if (channel == "Bpdpd0b") {
        std::vector<std::string> params = {
           "A1_dcu_Bpdpd0b_re",
	   "A1_dcu_Bpdpd0b_im"
        };
        channelParameters[channel] = params;
        AddParameter("A1_dcu_Bpdpd0b_re", -130., 150.);
	AddParameter("A1_dcu_Bpdpd0b_im", -130., 150.);
    } else if (channel == "Bpdpsd0b") {
        std::vector<std::string> params = {
	   "su3_A1_scu_Bpdpsd0b_re",
	   "su3_A1_scu_Bpdpsd0b_im",
        };
	 channelParameters[channel] = params;

	AddParameter("su3_A1_scu_Bpdpsd0b_re", -1., 1.);
	AddParameter("su3_A1_scu_Bpdpsd0b_im", -1., 1.);

    } else if (channel == "Bsdpsdms") {
        std::vector<std::string> params = {
            "su3_E1_scc_Bsdpsdms_re",
	    "su3_E1_scc_Bsdpsdms_im",
            "su3_A2_csc_Bsdmsdps_re",
	    "su3_A2_csc_Bsdmsdps_im",
            "su3_P1_GIM_sc_Bsdpsdms_re",
	    "su3_P1_GIM_sc_Bsdpsdms_im",
            "su3_P3_GIM_sc_Bsdpsdms_re",
	    "su3_P3_GIM_sc_Bsdpsdms_im"
        };
        channelParameters[channel] = params;
        for (const auto& param : params) {
            AddParameter(param, -1., 1.);
        }
    } else if (channel == "Bsdpdms") {
    std::vector<std::string> params = {
            "su3_E1_dcc_Bsdpdms_re",
	    "su3_E1_dcc_Bsdpdms_im",
            "su3_P1_GIM_dc_Bsdpdms_re",
	    "su3_P1_GIM_dc_Bsdpdms_im",

    };
    channelParameters[channel] = params;
    for (const auto& param : params) {
        AddParameter(param, -1., 1.);
    }
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
std::map<std::string, double> BqDqDqbar::DeclareParameters() {

    // Ensure channelParameters is populated
    for (const auto& channel : channelNames) {
        // Call DefineParameters once per channel
        if (channelParameters.find(channel) == channelParameters.end()) {
            BqDqDqbar::DefineParameters(channel);
        }

        // Add parameters to parameterValues
        auto it = channelParameters.find(channel);
        if (it != channelParameters.end()) {
            for (const auto& param : it->second) {
	      parameterValues[param] = 0; //initialize
	      cout << "parameter " << param << " added with value: " << parameterValues[param] << endl;
            }
        }
    }

    cout << "Declare Parameters called correctly" << endl;
    return parameterValues;
}

//--------------------------------------------------
Parameter BqDqDqbar::getPar(const std::string& baseName) const {
        // Look for the real and imaginary parts in the parameterValues map
    auto it_real = parameterValues.find(baseName + "_re");
    auto it_imag = parameterValues.find(baseName + "_im");

    if (it_real != parameterValues.end() && it_imag != parameterValues.end()) {
        // Construct a complex number from the real and imaginary parts
        TComplex value(it_real->second, it_imag->second);
        cout << "get par called correctly for par " << value << endl;
        return value;
    } else {
        throw std::runtime_error("Error: Real or imaginary part for parameter " + baseName + " not found.");
    }
}
//---------------------------------------------------------------------

double BqDqDqbar::getParameterValue(const std::string& paramName) const {
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



// Setter function: sets the value for a given parameter in the map
void BqDqDqbar::SetParameterValue(const std::string& paramName, double value) {
    parameterValues[paramName] = value;  // Insert or update the value for the given parameter name
}



 //compute decay amplitude for each channel
void BqDqDqbar::compute_decay_amplitudes(const std::string& channel, bool conjugate) {
   //cout << "computing decay amplitude for channel " << channel << endl;
   amplitudes[channel] = TComplex(0.0, 0.0);

    Parameter amp;

    // CKM elements with conjugation if necessary
    TComplex Vcd = conjugate ? ckm.getVdc() : ckm.getVcd();
    TComplex Vcs = conjugate ? ckm.getVsc() : ckm.getVcs();
    TComplex Vbc = conjugate ? ckm.getVcb() : ckm.getVbc();
    TComplex Vtd = conjugate ? ckm.getVdt() : ckm.getVtd();
    TComplex Vbt = conjugate ? ckm.getVtb() : ckm.getVbt();
    TComplex Vud = conjugate ? ckm.getVdu() : ckm.getVud();
    TComplex Vbu = conjugate ? ckm.getVub() : ckm.getVbu();
    TComplex Vts = conjugate ? ckm.getVst() : ckm.getVts();
    TComplex Vus = conjugate ? ckm.getVsu() : ckm.getVus();

    double A_real = 0.0;
    double phaseA = 0.0;
    double modA = 0.0;
    double A_imag = 0.0;
    double B_real = 0.0;
    double modB = 0.0;
    double phaseB = 0.0;
    double B_imag = 0.0;
    TComplex A(0.0, 0.0);
    TComplex B(0.0, 0.0);

 // Define the SU(3)-breaking for the E1 terms
TComplex E1_Bddpsdm = TComplex( getParameterValue("E1_dcc_Bddpdm_re") * (1 + getParameterValue("su3_E1_scc_Bddpsdm_re")),
    getParameterValue("E1_dcc_Bddpdm_re") * getParameterValue("su3_E1_scc_Bddpsdm_im")
);

TComplex E1_Bsdpsdms = TComplex(
							getParameterValue("E1_dcc_Bddpdm_re") * (1 + getParameterValue("su3_E1_scc_Bsdpsdms_re")),
							getParameterValue("E1_dcc_Bddpdm_re") * getParameterValue("su3_E1_scc_Bsdpsdms_im")
);

TComplex E1_Bsdpdms = TComplex(
						       getParameterValue("E1_dcc_Bddpdm_re") * (1 + getParameterValue("su3_E1_dcc_Bsdpdms_re")),
						       getParameterValue("E1_dcc_Bddpdm_re") * getParameterValue("su3_E1_dcc_Bsdpdms_im")
);


  TComplex P1_GIM_Bddpsdm = TComplex(
							     getParameterValue("P1_GIM_dc_Bddpdm_re") * (1 + getParameterValue("su3_P1_GIM_sc_Bddpsdm_re")),
							     getParameterValue("P1_GIM_dc_Bddpdm_im") * (1 +  getParameterValue("su3_P1_GIM_sc_Bddpsdm_im") )
);

  TComplex A1_Bpdpsd0b = TComplex(getParameterValue("A1_dcu_Bpdpd0b_re") * (1 + getParameterValue("su3_A1_scu_Bpdpsd0b_re")),
							  getParameterValue("A1_dcu_Bpdpd0b_im") * (1 +  getParameterValue("su3_A1_scu_Bpdpsd0b_im"))
);

  TComplex A2_Bsdpsdms = TComplex(getParameterValue("A2_cdc_Bddmdp_re") * (1+ getParameterValue("su3_A2_csc_Bsdmsdps_re")),
							  getParameterValue("A2_cdc_Bddmdp_im") * (1+ getParameterValue("su3_A2_csc_Bsdmsdps_im")) );


  TComplex P1_GIM_Bsdpsdms = TComplex(getParameterValue("P1_GIM_dc_Bddpdm_re") * (1+ getParameterValue("su3_P1_GIM_sc_Bsdpsdms_re")), getParameterValue("P1_GIM_dc_Bddpdm_im") * (1+ getParameterValue("su3_P1_GIM_sc_Bsdpsdms_im")));

  TComplex P3_GIM_Bsdpsdms = TComplex(getParameterValue("P3_GIM_dc_Bddmdp_re") * (1+ getParameterValue("su3_P3_GIM_sc_Bsdpsdms_re")),
							      getParameterValue("P3_GIM_dc_Bddmdp_im") * (1+ getParameterValue("su3_P3_GIM_sc_Bsdpsdms_im")) );


  TComplex P1_GIM_Bsdpdms = TComplex(getParameterValue("P1_GIM_dc_Bddpdm_re") * (1+getParameterValue("su3_P1_GIM_dc_Bsdpdms_re")),
							     getParameterValue("P1_GIM_dc_Bddpdm_im") * (1+getParameterValue("su3_P1_GIM_dc_Bsdpdms_im")));


    // Start computing the amplitude depending on the channel
    if (channel == "Bddpdm") {
        amp = Vcd * Vbc * (getPar("E1_dcc_Bddpdm") + getPar("A2_cdc_Bddmdp"))
	  - Vud * Vbu * (getPar("P1_GIM_dc_Bddpdm") + getPar("P3_GIM_dc_Bddmdp"));
        amplitudes[channel] = amp;
	A_real = getPar("E1_dcc_Bddpdm").Re() + getPar("A2_cdc_Bddmdp").Re();
	A_imag = getPar("A2_cdc_Bddmdp").Im();
	B_real = getPar("P1_GIM_dc_Bddpdm").Re() + getPar("P3_GIM_dc_Bddmdp").Re();
	B_imag = getPar("P1_GIM_dc_Bddpdm").Im() + getPar("P3_GIM_dc_Bddmdp").Im();
	A = getPar("E1_dcc_Bddpdm") + getPar("A2_cdc_Bddmdp");
	B = getPar("P1_GIM_dc_Bddpdm") + getPar("P3_GIM_dc_Bddmdp");
	phaseA = arg(A);
	phaseB = arg(B);
	modA = abs(A);
	modB = abs(B);





    } else if (channel == "Bddpsdm") {
      amp = Vcs * Vbc * (E1_Bddpsdm)
	- Vus * Vbu * (P1_GIM_Bddpsdm);
        amplitudes[channel] = amp;
	A_real = E1_Bddpsdm.Re();
	A_imag = E1_Bddpsdm.Im();
	B_real = P1_GIM_Bddpsdm.Re();
	B_imag = P1_GIM_Bddpsdm.Im();
	A = E1_Bddpsdm;
	B = P1_GIM_Bddpsdm;
	phaseA = arg(A);
	phaseB = arg(B);
	modA = abs(A);
	modB = abs(B);




    } else if (channel == "Bpdpd0b") {
      amp = Vcd * Vbc * getPar("E1_dcc_Bddpdm")
	+ Vud * Vbu * (getPar("A1_dcu_Bpdpd0b") - getPar("P1_GIM_dc_Bddpdm"));
        amplitudes[channel] = amp;
	A_real = getPar("E1_dcc_Bddpdm").Re();
	A_imag = getPar("E1_dcc_Bddpdm").Im();
	B_real = getPar("A1_dcu_Bpdpd0b").Re() - getPar("P1_GIM_dc_Bddpdm").Re();
	B_imag = getPar("A1_dcu_Bpdpd0b").Im() - getPar("P1_GIM_dc_Bddpdm").Im();
	A = getPar("E1_dcc_Bddpdm");
	B = getPar("A1_dcu_Bpdpd0b") - getPar("P1_GIM_dc_Bddpdm");
	phaseA = arg(A);
	phaseB = arg(B);
	modA = abs(A);
	modB = abs(B);




    } else if (channel == "Bpdpsd0b") {
      amp = Vcs * Vbc * (E1_Bddpsdm)
	  + Vus * Vbu * (A1_Bpdpsd0b - P1_GIM_Bddpsdm);
        amplitudes[channel] = amp;
	A_real = E1_Bddpsdm.Re();
	A_imag = E1_Bddpsdm.Im();
	B_real = A1_Bpdpsd0b.Re() - P1_GIM_Bddpsdm.Re();
	B_imag = A1_Bpdpsd0b.Im() - P1_GIM_Bddpsdm.Im();
	A = E1_Bddpsdm;
	B = A1_Bpdpsd0b - P1_GIM_Bddpsdm;
	phaseA = arg(A);
	phaseB = arg(B);
	modA = abs(A);
	modB = abs(B);




    } else if (channel == "Bsdpsdms") {
      amp = Vcs * Vbc * (E1_Bsdpsdms + A2_Bsdpsdms)
	  - Vus * Vbu * (P1_GIM_Bsdpsdms + P3_GIM_Bsdpsdms);
        amplitudes[channel] = amp;
	A_real = E1_Bsdpsdms.Re() + A2_Bsdpsdms.Re();
	A_imag = E1_Bsdpsdms.Im() + A2_Bsdpsdms.Im();
	B_real = P1_GIM_Bsdpsdms.Re() + P3_GIM_Bsdpsdms.Re();
	B_imag = P1_GIM_Bsdpsdms.Im() + P3_GIM_Bsdpsdms.Im();
	A = E1_Bsdpsdms + A2_Bsdpsdms;
	B = P1_GIM_Bsdpsdms + P3_GIM_Bsdpsdms;
	phaseA = arg(A);
	phaseB = arg(B);
	modA = abs(A);
	modB = abs(B);


    } else if (channel == "Bsdpdms") {
      amp = Vcd * Vbc * (E1_Bsdpdms)
	- Vud * Vbu * (P1_GIM_Bsdpdms);
        amplitudes[channel] = amp;

	A_real = E1_Bsdpdms.Re();
	A_imag = E1_Bsdpdms.Im();
	B_real = P1_GIM_Bsdpdms.Re();
	B_imag = P1_GIM_Bsdpdms.Im();
	A = E1_Bsdpdms;
	B = P1_GIM_Bsdpdms;
	phaseA = arg(A);
	phaseB = arg(B);
	modA = abs(A);
	modB = abs(B);



    }

    obs["A_" + channel + "_re"] = A_real;
    obs["A_" + channel + "_im"] = A_imag;
    obs["B_" + channel + "_re"] = B_real;
    obs["B_" + channel + "_im"] = B_imag;
    obs["|A|(" + channel + ")"] = modA;
    obs["A_phase (" + channel + ")"] = phaseA;
    obs["|B|(" + channel + ")"] = modB;
    obs["B_phase (" + channel +  ")"] = phaseB;


}

 //getter for amplitudes
Parameter BqDqDqbar::get_amplitude(const std::string& channel) {
  compute_decay_amplitudes(channel, false);  // false for normal
    return amplitudes[channel];
}


 //getter for conjugated amplitudes
Parameter BqDqDqbar::get_conjugate_amplitude(const std::string& channel) {
  compute_decay_amplitudes(channel, true);  // true for conjugate
    return amplitudes[channel];
}

// -----------------------------------------------------------
BqDqDqbar::~BqDqDqbar() {
  //delete histos;  // Clean up the dynamically allocated memory
};

// ---------------------------------------------------------
std::pair<std::string, std::pair<std::string, std::string>> BqDqDqbar::parseChannel(const std::string& channel) const {
    // First, identify the B meson part
    std::string bMeson;
    if (channel.substr(0, 2) == "Bp") {
        bMeson = "Bp";
    } else if (channel.substr(0, 2) == "Bd") {
        bMeson = "Bd";
    } else if (channel.substr(0, 2) == "Bs") {
        bMeson = "Bs";
    } else {
        throw std::runtime_error("Unknown B meson in channel: " + channel);
    }

    // Next, identify the final state mesons (everything after the first 2 letters)
    std::string remaining = channel.substr(2);
    std::string meson1, meson2;

    // Create a vector of meson names sorted by length in descending order
    std::vector<std::string> mesonNames;
    for (const auto& meson : mesonMasses) {
        mesonNames.push_back(meson.first);
    }
    std::sort(mesonNames.begin(), mesonNames.end(), [](const std::string& a, const std::string& b) {
        return a.length() > b.length();  // Sort by length in descending order
    });

    // Find the first final-state meson
    for (const auto& mesonName : mesonNames) {
        if (remaining.find(mesonName) == 0) {
            meson1 = mesonName;
            remaining = remaining.substr(mesonName.length());  // Remove meson1 from remaining
            break;
        }
    }

    // Find the second final-state meson
    for (const auto& mesonName : mesonNames) {
        if (remaining.find(mesonName) == 0) {
            meson2 = mesonName;
            break;
        }
    }

    if (meson1.empty() || meson2.empty()) {
        throw std::runtime_error("Unable to parse final-state mesons in channel: " + channel);
    }

    return {bMeson, {meson1, meson2}};
}





double BqDqDqbar::getBMesonLifetime(const std::string& bMeson) const {
    if (bMeson == "Bp") {
        return tau_Bp;
    } else if (bMeson == "Bd") {
        return tau_Bd;
    } else if (bMeson == "Bs") {
        return tau_Bs;
    } else {
        throw std::runtime_error("Unknown B meson: " + bMeson);
    }
}


double BqDqDqbar::getBMesonMass(const std::string& bMeson) const {
    if (bMeson == "Bp") {
        return m_Bp;
    } else if (bMeson == "Bd") {
        return m_Bd;
    } else if (bMeson == "Bs") {
        return m_Bs;
    } else {
        throw std::runtime_error("Unknown B meson: " + bMeson);
    }
}



double BqDqDqbar::CalculateBR(Parameter amplitude, const string& channel) const {
   double BR;
   // Parse the channel to get the decaying B meson and final-state mesons
    auto parsed = parseChannel(channel);
    std::string bMeson = parsed.first;
    std::string meson1 = parsed.second.first;
    std::string meson2 = parsed.second.second;

    // Get the mass of the decaying B meson
    double m_B = getBMesonMass(bMeson);

    // Get the masses of the final-state mesons
    double m1 = mesonMasses.at(meson1);
    double m2 = mesonMasses.at(meson2);


    // Calculate the momentum magnitude |p| of the final-state mesons
    double p = sqrt((m_B * m_B - (m1 + m2) * (m1 + m2)) * (m_B * m_B - (m1 - m2) * (m1 - m2))) / (2.0 * m_B);

    // The decay width (Γ) is proportional to |p| and the square of the amplitude
    double decay_width = (std::pow(G_F, 2) / 2) * (p / (std::pow(m_B, 2) * 16 * M_PI * h_t)) * (std::norm(amplitude));


    // Get the lifetime and calculate total width
    double lifetime = getBMesonLifetime(bMeson);
    BR = decay_width*lifetime;
    // Check if the final-state mesons are identical and apply the 1/2 factor
    if (meson1 == meson2) {
        BR *= 0.5;
    }

    return BR;

 }

///-----------------------------------------------------------------------
 // Function to calculate A_CP asymmetry
double BqDqDqbar::CalculateAcp(const Parameter& amplitude, const Parameter& conjugate_amplitude) const {
    // Calculate the squared norms of the amplitudes
    double A2 = std::norm(amplitude);
    double Abar2 = std::norm(conjugate_amplitude);

    // Compute A_CP
    double Acp = (Abar2 - A2) / (A2 + Abar2);

    return Acp;
}
//---------------------------------------------------------------------
double BqDqDqbar::CalculateC(const Parameter& amplitude, const Parameter& conjugate_amplitude, const std::string& channel) {
    // Get q/p based on whether the channel is Bd or Bs (use parseChannel if necessary)
    auto parsed = parseChannel(channel);
    std::string bMeson = parsed.first;

    TComplex q_p = (bMeson == "Bd") ? ckm.get_q_p_Bd() : ckm.get_q_p_Bs();

    // Calculate the ratio lambda = eta * (q/p) * (A_cp / A_conj)
    double eta = cpEigenvalue[channel];  // Assume cpEigenvalue map has eta for each channel
    TComplex lambda = eta * q_p * (conjugate_amplitude / amplitude);

    // Calculate C observable: C = (1 - |lambda|^2) / (1 + |lambda|^2)
    double mod_lambda_squared = std::norm(lambda);
    double C = (1.0 - mod_lambda_squared) / (1.0 + mod_lambda_squared);

    return C;
}

 double BqDqDqbar::CalculateS(const Parameter& amplitude, const Parameter& conjugate_amplitude, const std::string& channel) {
    // Get q/p based on whether the channel is Bd or Bs
    auto parsed = parseChannel(channel);
    std::string bMeson = parsed.first;

    TComplex q_p = (bMeson == "Bd") ? ckm.get_q_p_Bd() : ckm.get_q_p_Bs();

    // Calculate the ratio lambda = eta * (q/p) * (A_cp / A_conj)
    double eta = cpEigenvalue[channel];  // Assume cpEigenvalue map has eta for each channel
    TComplex lambda = eta * q_p * (conjugate_amplitude / amplitude);

    // Calculate S observable: S = 2 Im(lambda) / (1 + |lambda|^2)
    double mod_lambda_squared = std::norm(lambda);
    double S = 2.0 * std::imag(lambda) / (1.0 + mod_lambda_squared);

    return S;
}




std::pair<double, double> BqDqDqbar::CalculatePhiAndLambda(const Parameter& amplitude, const Parameter& conjugate_amplitude, const std::string& channel) {
    // Get q/p based on the B meson (Bs in this case)
    TComplex q_p = ckm.get_q_p_Bs();

    // Compute lambda = (q/p) * (A_cp / A_conj)
    TComplex lambda = q_p * (conjugate_amplitude / amplitude);

    // Compute |lambda|
    double mod_lambda = std::abs(lambda);

    // Compute phi_s = -arg(lambda)
    double phi_s = -std::arg(lambda);

    return {phi_s, mod_lambda};
}



std::pair<std::vector<std::string>, std::string> BqDqDqbar::extractChannelFromCorrKey(const std::string& corr_key) {
    std::vector<std::string> channels;
    std::string experiment;

    // First, check if the key starts with "CS_" for C and S observables
    if (corr_key.find("CS_") == 0) {
        // Key format: "CS_channel_exp"
        size_t underscore = corr_key.find("_", 3);  // Find underscore after "CS_"

        // Extract channel and experiment
        std::string channel = corr_key.substr(3, underscore - 3);  // Extract channel name
        experiment = corr_key.substr(underscore + 1);              // Extract experiment name

        // Store the single channel in the vector
        channels.push_back(channel);

    } else if (corr_key.find("ACP_") == 0) {
        // Key format: "ACP_channel1_channel2_exp"
        size_t first_underscore = corr_key.find("_", 4);  // Find underscore after "ACP_"
        size_t second_underscore = corr_key.find("_", first_underscore + 1);

        // Extract the two channels and the experiment
        std::string channel1 = corr_key.substr(4, first_underscore - 4);  // Extract first channel
        std::string channel2 = corr_key.substr(first_underscore + 1, second_underscore - first_underscore - 1); // Extract second channel
        experiment = corr_key.substr(second_underscore + 1);  // Extract experiment name

        // Store both channels in the vector
        channels.push_back(channel1);
        channels.push_back(channel2);

    } else if (corr_key.find("phi_lamda_") == 0) {
        // Key format: "phi_lamda_channel_exp"
        size_t underscore = corr_key.find("_", 10);  // Find underscore after "phi_lamda_"

        // Extract channel and experiment
        std::string channel = corr_key.substr(10, underscore - 10);  // Extract channel name
        experiment = corr_key.substr(underscore + 1);                // Extract experiment name

        // Store the single channel in the vector
        channels.push_back(channel);
    } else {
        throw std::runtime_error("Unknown key format: " + corr_key);
    }

    return {channels, experiment};
}


//------------------------------------------------------------------------
double BqDqDqbar::Calculate_UncorrelatedObservables(const std::map<std::string, std::pair<Parameter, Parameter>>& amplitude_map) {
    double ll_uncorr = 0.0;  // Initialize the log-likelihood contribution


    // ACP for Bpdpd0b
    {
        const auto& amp_pair_Bpdpd0b = amplitude_map.at("Bpdpd0b");
        double acp_predicted = CalculateAcp(amp_pair_Bpdpd0b.first, amp_pair_Bpdpd0b.second);
        obs["ACP_Bpdpd0b"] = acp_predicted;

        // Get measurement for ACP
        try {
            double acp_observed = meas.at("ACPBpdpd0b").getMean();
            double acp_uncertainty = meas.at("ACPBpdpd0b").getSigma();
            double diff_acp = acp_predicted - acp_observed;
            ll_uncorr += -0.5 * (diff_acp * diff_acp / (acp_uncertainty * acp_uncertainty));
        } catch (const std::out_of_range&) {
            std::cerr << "Error: ACP for Bpdpd0b not found in meas map." << std::endl;
        }
    }

    // C and S for Bddpdm
    {
        const auto& amp_pair_Bddpdm = amplitude_map.at("Bddpdm");
        double c_predicted = CalculateC(amp_pair_Bddpdm.first, amp_pair_Bddpdm.second, "Bddpdm");
        double s_predicted = CalculateS(amp_pair_Bddpdm.first, amp_pair_Bddpdm.second, "Bddpdm");
        obs["C_Bddpdm"] = c_predicted;
        obs["S_Bddpdm"] = s_predicted;

        // Get measurements for C and S
        try {
            double c_observed = meas.at("CBddpdm").getMean();
            double c_uncertainty = meas.at("CBddpdm").getSigma();
            double diff_c = c_predicted - c_observed;
            ll_uncorr += -0.5 * (diff_c * diff_c / (c_uncertainty * c_uncertainty));

            double s_observed = meas.at("SBddpdm").getMean();
            double s_uncertainty = meas.at("SBddpdm").getSigma();
            double diff_s = s_predicted - s_observed;
            ll_uncorr += -0.5 * (diff_s * diff_s / (s_uncertainty * s_uncertainty));
        } catch (const std::out_of_range&) {
            std::cerr << "Error: C or S for Bddpdm not found in meas map." << std::endl;
        }
    }

    //mod(lambda) contribution to predict phi_s
    {
        const auto& amp_pair_Bsdpsdms = amplitude_map.at("Bsdpsdms");
        auto phi_lambda_result = CalculatePhiAndLambda(amp_pair_Bsdpsdms.first, amp_pair_Bsdpsdms.second, "Bsdpsdms");
	double lamda_predicted = phi_lambda_result.second;
	double phi_predicted = phi_lambda_result.first;
	double beta_s = arg(-ckm.getVts() * ckm.getVbt() / (ckm.getVcs() * ckm.getVbc()));
	double phi_th = -2*beta_s;
	double delta_phi = phi_predicted - phi_th;
	obs["delta_phi"] = delta_phi;

	// Store phi_s and |lambda| in the obs map
        obs["phi_Bsdpsdms"] = phi_lambda_result.first;   // phi_s for Bsdpsdms
        obs["lamda_Bsdpsdms"] = phi_lambda_result.second;  // |lambda| for Bsdpsdms


    try {
      double lamda_observed = meas.at("L_Bsdpsdms").getMean();
      double lamda_uncertainty = meas.at("L_Bsdpsdms").getSigma();
      double diff_L = lamda_predicted - lamda_observed;
      ll_uncorr += -0.5 * (diff_L * diff_L / (lamda_uncertainty * lamda_uncertainty));
    } catch (const std::out_of_range&) {
            std::cerr << "Error: mod(lamda) not found in meas map." << std::endl;
        }
    }


    return ll_uncorr;
}


 //----------------------------------------------------------------------------------

double BqDqDqbar::Calculate_CorrelatedObservables(const std::map<std::string, std::pair<Parameter, Parameter>>& amplitude_map) {
    double ll_corr = 0.0;  // Initialize the log-likelihood contribution
    TVectorD corr(2);  // For correlated observables (e.g., C and S)

    // **Correlated observables for Bddpdm**
    {
        const auto& amp_pair_Bddpdm = amplitude_map.at("Bddpdm");
        corr(0) = CalculateC(amp_pair_Bddpdm.first, amp_pair_Bddpdm.second, "Bddpdm");
        corr(1) = CalculateS(amp_pair_Bddpdm.first, amp_pair_Bddpdm.second, "Bddpdm");
        ll_corr += corrmeas.at("CS_Bddpdm_BaBar2008").logweight(corr);  // BaBar:2008xnt contribution
        ll_corr += corrmeas.at("CS_Bddpdm_LHCb2024").logweight(corr);   // LHCb:2024gkk contribution
    }

    // **Correlated observables for Bpdpd0b and Bpdpsd0b**
    {
        const auto& amp_pair_Bpdpd0b = amplitude_map.at("Bpdpd0b");
        const auto& amp_pair_Bpdpsd0b = amplitude_map.at("Bpdpsd0b");


        // Calculate and store ACP for Bpdpd0b
        double acp_Bpdpd0b = CalculateAcp(amp_pair_Bpdpd0b.first, amp_pair_Bpdpd0b.second);
        corr(0) = acp_Bpdpd0b;

        // Calculate and store ACP for Bpdpsd0b
        double acp_Bpdpsd0b = CalculateAcp(amp_pair_Bpdpsd0b.first, amp_pair_Bpdpsd0b.second);
        obs["ACP_Bpdpsd0b"] = acp_Bpdpsd0b;  // Store in obs map for histogram
        corr(1) = acp_Bpdpsd0b;

        ll_corr += corrmeas.at("ACP_Bpdpd0b_Bpdpsd0b_LHCb2023").logweight(corr);
    }

    return ll_corr;
}


//------------------------------------------------------------

double BqDqDqbar::LogLikelihood(const std::vector<double>& parameters) {
    static int iteration_counter = 0;  // Static to persist across calls
    ++iteration_counter;  // Increment the counter on every call
    if (iteration_counter % 50000 == 0) {
        cout << "Iteration " << iteration_counter << " processing... " << endl;
    }

    obs.clear();  // Clear obs map for each iteration
    double ll = 0.0;  // Log-likelihood accumulator

    int expectedSize = 0;
    for (const auto& channel : channelNames) {
      auto it = channelParameters.find(channel);
      if (it != channelParameters.end()) {
	expectedSize += it->second.size() ;  // Each complex parameter has 2 values (real, imaginary)
      }
    }
    //std::cout << "Expected parameters size: " << expectedSize << std::endl;
    if (parameters.size() != expectedSize + 4) {
      std::cerr << "Error: parameters.size() = " << parameters.size()
              << ", but expectedSize = " << expectedSize + 4 << std::endl;
      return -1e30;
    }

    // Ensure that the parameters are valid before using them
    if (parameters.empty()) {
        std::cerr << "Error: Empty parameters vector!" << std::endl;
        return -1e30;
    }

    // Unpack CKM parameters from the last 4 entries in the parameters vector
    std::vector<double> ckmParams(parameters.end() - 4, parameters.end());
     //compute ckm elements
    ckm.computeCKM(ckmParams[0], ckmParams[1], ckmParams[2], ckmParams[3], true);


    // Unpack the parameters and store them in the parametersValue map
    std::map<std::string, std::pair<Parameter, Parameter>> amplitude_map;
    amplitude_map.clear();  // Ensure the map is cleared for every iteration

       // Unpack the parameters and store them in the parametersValue map
   int index = 0; // This index tracks our position in the 'parameters' vector
   for (const auto& channel : channelNames) { // Loop through each channel name
     auto it = channelParameters.find(channel); // Find parameters for the current channel
     if (it != channelParameters.end()) { // Ensure the channel is found in channelParameters
        // Loop over parameters of this channel, in steps of 2 (since params are stored as real, imaginary pairs)
       for (size_t j = 0; j < it->second.size(); j += 2) {
            // Ensure that we don't go out of bounds when accessing the parameters vector
	 if (index + 1 >= parameters.size() - 4 ) {
	   break;
	 }

            // Extract real and imaginary parts from the 'parameters' vector at the current index
         double realPart = parameters[index];       // Extract the real part
         double imagPart = parameters[index + 1];   // Extract the imaginary part

            // Get the corresponding parameter names for the real and imaginary parts from the channel
         const std::string& realParamName = it->second[j];     // Name for the real part
         const std::string& imagParamName = it->second[j + 1]; // Name for the imaginary part

            // Use the setter method to store these values into 'parameterValues' map
         SetParameterValue(realParamName, realPart);  // Store real part
         SetParameterValue(imagParamName, imagPart);  // Store imaginary part
	 //cout << "processing channel = " << channel << endl;

            // Increment the index by 2, because we're handling two values (real, imaginary) at a time
         index += 2;
       }
     } else {
        // If the channel is not found in 'channelParameters', print an error message and return a failure code
       std::cerr << "Channel " << channel << " not found in channelParameters." << std::endl;
       return -1e30; // Return an error code
     }
   }

   // Loop over channelNames

    for (const std::string& channel : channelNames) {
     // Compute and store amplitudes for each channel
      amplitude_map[channel] = std::make_pair(
                get_amplitude(channel),
                get_conjugate_amplitude(channel)
            );
    //Step 1: Use amplitudes to calculate BR and log-likelihood contributions

        const auto& amp_pair = amplitude_map[channel];  // Get precomputed amplitude and conjugate amplitude

        if (std::isnan(amp_pair.first.Re()) || std::isnan(amp_pair.first.Im()) ||
            std::isinf(amp_pair.first.Re()) || std::isinf(amp_pair.first.Im()) ||
            std::isnan(amp_pair.second.Re()) || std::isnan(amp_pair.second.Im()) ||
            std::isinf(amp_pair.second.Re()) || std::isinf(amp_pair.second.Im())) {
            std::cerr << "Invalid amplitude (NaN or Inf) for channel: " << channel << std::endl;
            return -100;
        }


    //construct histograms for effective parameters
    Parameter E1 = getPar("E1_dcc_Bddpdm");
	double E1_re = getPar("E1_dcc_Bddpdm").Re();

    double modE1 = TComplex::Abs(E1);
    double phaseE1 = TComplex::Arg(E1);

    Parameter A2 = getPar("A2_cdc_Bddmdp");
	double A2_re = getPar("A2_cdc_Bddmdp").Re();
	double A2_im = getPar("A2_cdc_Bddmdp").Im();

    double modA2 = TComplex::Abs(A2);
    double phaseA2 = TComplex::Arg(A2);


    Parameter A1 = getPar("A1_dcu_Bpdpd0b");
	double A1_re = getPar("A1_dcu_Bpdpd0b").Re();
	double A1_im = getPar("A1_dcu_Bpdpd0b").Im();


    double modA1 = TComplex::Abs(A1);
    double phaseA1 = TComplex::Arg(A1);


    Parameter P1_GIM = getPar("P1_GIM_dc_Bddpdm");
	double P1_GIM_re = getPar("P1_GIM_dc_Bddpdm").Re();
	double P1_GIM_im = getPar("P1_GIM_dc_Bddpdm").Im();


    double modP1_GIM = TComplex::Abs(P1_GIM);
    double phaseP1_GIM = TComplex::Arg(P1_GIM);


    Parameter P3_GIM = getPar("P3_GIM_dc_Bddmdp");
	double P3_GIM_re = getPar("P3_GIM_dc_Bddmdp").Re();
	double P3_GIM_im = getPar("P3_GIM_dc_Bddmdp").Im();

    double modP3_GIM = TComplex::Abs(P3_GIM);
    double phaseP3_GIM = TComplex::Arg(P3_GIM);








	obs["|E1|"] = modE1;
	obs["|A2|"] = modA2;
	obs["phase A2"] = phaseA2;
	obs["|A1|"] = modA1;
	obs["phase A1"] = phaseA1;
	obs["|P1_GIM|"] = modP1_GIM;
	obs["phase P1_GIM"] = phaseP1_GIM;
	obs["|P3_GIM|"] = modP3_GIM;
	obs["phase P3_GIM"] = phaseP3_GIM;
	obs["E1"] = E1_re;
	obs["A2_re"] = A2_re;
	obs["A2_im"] = A2_im;
	obs["A1_re"] = A1_re;
	obs["A1_im"] = A1_im;
	obs["P1_GIM_re"] = P1_GIM_re;
	obs["P1_GIM_im"] = P1_GIM_im;
	obs["P3_GIM_re"] = P3_GIM_re;
	obs["P3_GIM_im"] = P3_GIM_im;




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

    //Step 2: Calculate the log-likelihood contribution for uncorrelated observables
    ll += Calculate_UncorrelatedObservables(amplitude_map);

    //Step 3: Calculate the log-likelihood contribution for correlated observables
    ll += Calculate_CorrelatedObservables(amplitude_map);


    return ll;
}


//---------------------------------------------------------


void BqDqDqbar::MCMCUserIterationInterface() {
    // Loop over all MCMC chains
  const unsigned int log_interval = 10000;  // Log every 1000 iterations
  std::vector<double> pars;

    for (unsigned int i = 0; i < fMCMCNChains; ++i) {
        pars = fMCMCStates.at(i).parameters;
        LogLikelihood(pars);

        // Directly fill histograms
        histos.fillh1d();
        histos.fillh2d();
    }
}



void BqDqDqbar::SaveHistograms(const std::string& filename) {
    TFile file(filename.c_str(), "RECREATE");

    // Write all histograms to the file
    for (const auto& histPair : histos.h1d) {
        histPair.second->Write();
    }

    file.Close();
    std::cout << "Histograms saved to " << filename << std::endl;
}

// ---------------------------------------------------------

void BqDqDqbar::PrintObservablePulls(const std::string& filename) {
    std::ofstream outfile(filename);

    if (!outfile.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return;
    }

    // Header for the output file
    outfile << "Observable\tMeasurement\tMean\tSigma\tPull\n";

    // Loop through each 1D histogram for observables
    for (const auto& histPair : histos.h1d) {
        const std::string& obs_name = histPair.first;
        TH1D* hist = histPair.second;

        // Skip specific observables that don't require a pull calculation
        if (obs_name == "delta_phi" || obs_name == "phi_Bsdpsdms") {
            continue;
        }

        // Calculate the mean value from the histogram
        double obs_mean = hist->GetMean();

        // Manually map each observable to its corresponding measurement in meas
        double obs_measurement = 0.0;
        double sigma_measurement = 0.0;
        bool found_measurement = false;

        if (obs_name == "BR_Bddpdm" && meas.find("BRBddpdm") != meas.end()) {
          obs_measurement = meas.at("BRBddpdm").getMean();
          sigma_measurement = meas.at("BRBddpdm").getSigma();
          found_measurement = true;
        } else if (obs_name == "BR_Bddpsdm" && meas.find("BRBddpsdm") != meas.end()) {
          obs_measurement = meas.at("BRBddpsdm").getMean();
          sigma_measurement = meas.at("BRBddpsdm").getSigma();
          found_measurement = true;
        } else if (obs_name == "BR_Bpdpd0b" && meas.find("BRBpdpd0b") != meas.end()) {
	  obs_measurement = meas.at("BRBpdpd0b").getMean();
	  sigma_measurement = meas.at("BRBpdpd0b").getSigma();
	  found_measurement = true;
	} else if (obs_name == "BR_Bpdpsd0b" && meas.find("BRBpdpsd0b") != meas.end()) {
	  obs_measurement = meas.at("BRBpdpsd0b").getMean();
	  sigma_measurement = meas.at("BRBpdpsd0b").getSigma();
	  found_measurement = true;
	} else if (obs_name == "BR_Bsdpsdms" && meas.find("BRBsdpsdms") != meas.end()) {
	  obs_measurement = meas.at("BRBsdpsdms").getMean();
	  sigma_measurement = meas.at("BRBsdpsdms").getSigma();
	  found_measurement = true;
	} else if (obs_name == "BR_Bsdpdms" && meas.find("BRBsdpdms") != meas.end()) {
	  obs_measurement = meas.at("BRBsdpdms").getMean();
	  sigma_measurement = meas.at("BRBsdpdms").getSigma();
	  found_measurement = true;
	} else if (obs_name == "C_Bddpdm" && newmeas.find("CBddpdm") != meas.end()) {
	  obs_measurement = newmeas.at("CBddpdm").getMean();
          sigma_measurement = newmeas.at("CBddpdm").getSigma();
          found_measurement = true;
        } else if (obs_name == "S_Bddpdm" && newmeas.find("SBddpdm") != meas.end()) {
	  obs_measurement = newmeas.at("SBddpdm").getMean();
          sigma_measurement = newmeas.at("SBddpdm").getSigma();
          found_measurement = true;
        } else if (obs_name == "ACP_Bpdpd0b" && newmeas.find("ACPBpdpd0b") != meas.end()) {
	  obs_measurement = newmeas.at("ACPBpdpd0b").getMean();
          sigma_measurement = newmeas.at("ACPBpdpd0b").getSigma();
          found_measurement = true;
	} else if (obs_name == "ACP_Bpdpsd0b") {
	  obs_measurement = newmeas.at("ACPBpdpsd0b").getMean();
	  sigma_measurement = newmeas.at("ACPBpdpsd0b").getSigma();
	  found_measurement = true;
	} else if (obs_name == "lamda_Bsdpsdms") {
	  obs_measurement = meas.at("L_Bsdpsdms").getMean();
	  sigma_measurement = meas.at("L_Bsdpsdms").getSigma();
	  found_measurement = true;
	}

        if (found_measurement) {
            // Calculate the pull value
            // Create a Pull instance with required arguments
            Pull pullCalculator(
				*hist,                   // Reference to the histogram
				hist->GetNbinsX(),       // Number of bins along X
				500,                     // Number of bins along Y (arbitrary for now)
				0.0,                     // X-axis lower bound, triggers auto setting in makeCompatPlot
				0.0,                     // X-axis upper bound, triggers auto setting in makeCompatPlot
				0.0,                     // Y-axis lower bound, triggers auto setting in makeCompatPlot
				0.0,                     // Y-axis upper bound, triggers auto setting in makeCompatPlot
				true                     // Whether to use low stats or not
				);


            double pull = pullCalculator.calcPull(obs_measurement, sigma_measurement, true);
            outfile << obs_name << "\t"
                    << obs_measurement << "\t"
	            << sigma_measurement << "\t"
                    << obs_mean << "\t"
                    << pull << "\n";
        } else {
	  continue;
        }
    }

    outfile.close();
    std::cout << "Pull values saved to " << filename << std::endl;

}
