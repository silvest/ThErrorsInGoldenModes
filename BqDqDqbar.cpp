#include "BqDqDqbar.h"
#include "dato.h"
#include "PDGAverage.h"
#include "CKM.h"

#include <BAT/BCMath.h>
#include <BAT/BCGaussianPrior.h>
#include <TComplex.h>
#include <TRandom.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <complex>
#include <vector>
#include <algorithm>  // for std::find_if
#include <cctype>     // for std::isalpha
#include <stdexcept>  // for std::runtime_error
#include <TRandom3.h>
#include <cmath>
#include <TMatrixD.h>
#include <TDecompChol.h>
#include <iostream>
#include <fstream>
using namespace std;

// ---------------------------------------------------------

CKMParameters ckm;

BqDqDqbar::BqDqDqbar() : BCModel() //, histos(obs)
{
  cout << "constructor for BqDqDqbar called: inserting the experimental data" << endl;
    std::vector<dato> data;
    PDGAverage pdgaverage;
    //std::map<std::string, dato> meas;  // Map for measurements

    
    DeclareParameters();  // Ensure parameters are defined
     ckm.sampleCKMParameters();  

    

	
  

    /*for (const auto& channel : channelNames) {
      auto it = channelParameters.find(channel);
     if (it != channelParameters.end()) {
     for (size_t j = 0; j < it->second.size(); j += 2) {
               // Create histograms for real, imaginary, norm, and phase parts
                histos->createH1D(channel + "_re", 100, -10, 10);    // Real part histogram
                histos->createH1D(channel + "_im", 100, -10, 10);    // Imaginary part histogram
                histos->createH1D(channel + "_norm", 100, 0, 10);    // Norm histogram
                histos->createH1D(channel + "_phase", 100, -M_PI, M_PI);  // Phase histogram
            }
            }
        }
    */
   

    //measurments

    //BRBddpdm
    data.push_back(dato(2.12e-4, 0.16e-4, 0.18e-4)); //Belle:2012mef
    data.push_back(dato(1.97e-4, 0.20e-4, 0.20e-4)); //Belle:2007ebz
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

    //Bpdpd0bVSBpdpsd0b
    //meas.insert(pair<string, dato>("Bpdpd0bVSBpdpsd0b",dato(7.25e-2, 0.09e-2, 0.09e-2 ))); //LHCb:2023wbb


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
    data.push_back(dato(-0.43, 0.16, 0.05)); //Belle:2012mef
    data.push_back(dato(-0.91, 0.23, 0.06)); //Belle:2007ebz
    

    pdgaverage.setData(data);
    pdgaverage.setName("CBddpdm");
    pdgaverage.CalculateAverage();
    
    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    //SBddpdm
    data.push_back(dato(-0.99, 0.175, 0.08)); //Belle:2012mef
    data.push_back(dato(-1.13, 0.37, 0.09)); //Belle:2007ebz
    

    pdgaverage.setData(data);
    pdgaverage.setName("SBddpdm");
    pdgaverage.CalculateAverage();
    
    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    //ACPBpdpd0b
    data.push_back(dato(0e-2, 8e-2, 22-2)); //Belle:2008doh
    data.push_back(dato(-13e-2, 1e-24, 2e-2));  //BaBar:2006uih
    

    pdgaverage.setData(data);
    pdgaverage.setName("ACPBpdpd0b");
    pdgaverage.CalculateAverage();
    
    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();



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
    //CorrStat.Clear();
    //CorrSyst.Clear();
   

    //ACPBpdpd0b and ACPBpdpsd0b from LHCb:2023wbb
    dato Bpdpd0b(2.5e-2, 1.0e-2, 0.4e-2, 0.3e-2);
    dato Bpdpsd0b(0.5e-2, 0.2e-2, 0.5e-2, 0.3e-2);

    CorrData.push_back(dato(Bpdpd0b.getMean(), Bpdpd0b.getSigma()));  // Bpdpd0b
    CorrData.push_back(dato(Bpdpsd0b.getMean(), Bpdpsd0b.getSigma())); // Bpdpsd0b

    TMatrixDSym Corr(2); //total correlation

    
// Populate the correlation matrix
    Corr(0, 0) = 1.0;       // Variance for C
    Corr(1, 1) = 1.0;       // Variance for S
    Corr(0, 1) = 0.386;      // Correlation between C and S
    Corr(1, 0) = 0.386;      // Symmetric part (same as Corr(0, 1))


// Insert correlated data for LHCb:2016inx into corrmeas
    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("ACP_Bpdpd0b_Bpdpsd0b_LHCb2023",CorrelatedGaussianObservables(CorrData, Corr)));

    CorrData.clear();
    //Corr.Clear();

    //Bsdpsdms : modlambda phis from LHCb:2024gkk
    CorrData.push_back(dato(-0.055, 0.090, 0.021));  // phi_s 
    CorrData.push_back(dato(1.054, 0.099, 0.020)); // |lamda| observable 


// Populate the correlation matrix
    CorrStat(0, 0) = 1.;       // Variance for phi
    CorrStat(1, 1) = 1.;       // Variance for |lamda|
    CorrStat(0, 1) = 0.005;      // Correlation between phi and lamda
    CorrStat(1, 0) = 0.005;      // Symmetric part (same as Corr(0, 1))

    // Populate the correlation matrix
    CorrSyst(0, 0) = 1.;       //
    CorrSyst(1, 1) = 1.;  
    CorrSyst(0, 1) = 0.;      // Correlation between C and S
    CorrSyst(1, 0) = 0.;      // Symmetric part (same as Corr(0, 1))

// Insert correlated data for LHCb:2024gkk into corrmeas
    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("phi_lamda_Bsdpsdms_LHCb2024",CorrelatedGaussianObservables(CorrData, CorrStat, CorrSyst)));

    CorrData.clear();
    //CorrStat.Clear();
    //CorrSyst.Clear();
    std::cout << "All meas inserted" << endl;
    
}
    
    
 
    //---------------------------------------------------------------


    
  //method to define the parameters needed to calculate each decay amplitude, uses BCModel method AddParameter
void BqDqDqbar::DefineParameters(const string& channel) {
   double limit1 = 0.;
   //double limit2 = 1.5; 

  if (channel == "Bddpdm") {
        std::vector<std::string> params = {
            "E1_dcc_Bddpdm_re",
	    "E1_dcc_Bddpdm_im",
            "A2_cdc_Bddmdp_re",
	    "A2_cdc_Bddmdp_im",
            "P1_GIM_dc_Bddpdm_re",
	    "P1_GIM_dc_Bddpdm_im",
            "P3_GIM_dc_Bddmdp_re",
	    "P3_GIM_dc_Bddmdp_im",
            "P1_dc_Bddpdm_re",
	    "P1_dc_Bddpdm_im",
            "P3_dc_Bddmdp_re",
	    "P3_dc_Bddmdp_im"
        };
        channelParameters[channel] = params;
        for (const auto& param : params) {
            AddParameter(param, limit1, 50.);
        }
    } 
  else if (channel == "Bddpsdm") {
        std::vector<std::string> params = {
            "E1_scc_Bddpsdm_re",
	    "E1_scc_Bddpsdm_im",
            "P1_GIM_sc_Bddpsdm_re",
	    "P1_GIM_sc_Bddpsdm_im",
            "P1_sc_Bddpsdm_re",
	    "P1_sc_Bddpsdm_im"
	};
        channelParameters[channel] = params;
        for (const auto& param : params) {
            AddParameter(param, limit1, 50.);
        }
    } else if (channel == "Bpdpd0b") {
        std::vector<std::string> params = {
            "E1_dcc_Bpdpd0b_re",
	    "E1_dcc_Bpdpd0b_im",
            "P1_dc_Bpdpd0b_re",
	    "P1_dc_Bpdpd0b_im",
            "A1_dcu_Bpdpd0b_re",
	    "A1_dcu_Bpdpd0b_im",
            "P1_GIM_dc_Bpdpd0b_re",
	    "P1_GIM_dc_Bpdpd0b_im"
        };
        channelParameters[channel] = params;
        for (const auto& param : params) {
            AddParameter(param, limit1, 50.);
        }
    } else if (channel == "Bpdpsd0b") {
        std::vector<std::string> params = {
            "E1_scc_Bpdpsd0b_re",
	    "E1_scc_Bpdpsd0b_im",
            "P1_sc_Bpdpsd0b_re",
	    "P1_sc_Bpdpsd0b_im",
            "A1_scu_Bpdpsd0b_re",
	    "A1_scu_Bpdpsd0b_im",
            "P1_GIM_sc_Bpdpsd0b_re",
	    "P1_GIM_sc_Bpdpsd0b_im"
        };
        channelParameters[channel] = params;
        for (const auto& param : params) {
            AddParameter(param, limit1, 50.);
        }
    } else if (channel == "Bsdpsdms") {
        std::vector<std::string> params = {
            "E1_scc_Bsdpsdms_re",
	    "E1_scc_Bsdpsdms_im",
            "A2_csc_Bsdmsdps_re",
	    "A2_csc_Bsdmsdps_im",
            "P1_sc_Bsdpsdms_re",
	    "P1_sc_Bsdpsdms_im",
            "P3_sc_Bsdpsdms_re",
	    "P3_sc_Bsdpsdms_im",
            "P1_GIM_sc_Bsdpsdms_re",
	    "P1_GIM_sc_Bsdpsdms_im",
            "P3_GIM_sc_Bsdpsdms_re",
	    "P3_GIM_sc_Bsdpsdms_im"
	    
        };
        channelParameters[channel] = params;
        for (const auto& param : params) {
            AddParameter(param, limit1, 50.);
        }
    } else if (channel == "Bsdpdms") {
    std::vector<std::string> params = {
            "E1_dcc_Bsdpdms_re",
	    "E1_dcc_Bsdpdms_im",
            "P1_dc_Bsdpdms_re",
	    "P1_dc_Bsdpdms_im",
            "P1_GIM_dc_Bsdpdms_re",
	    "P1_GIM_dc_Bsdpdms_im"
    };
    channelParameters[channel] = params;
    for (const auto& param : params) {
        AddParameter(param, limit1, 10.);
    }
    } else {
        // Handle the case where the channel is not recognized
        std::cerr << "Error: Unrecognized channel \"" << channel << "\" in DefineParameters." << std::endl;
        throw std::runtime_error("Unrecognized channel: " + channel);
	} 
  std::cout << "Number of parameters: " << GetNParameters() << std::endl;

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
        } else {
            std::cerr << "Channel " << channel << " not found in channelParameters." << std::endl;
        }
    }
    cout << "Declare Parameters called correctly" << endl;
    return parameterValues;
}

Parameter BqDqDqbar::getPar(const std::string& baseName) const {
    // Look for the real and imaginary parts in the parameterValues map
    auto it_real = parameterValues.find(baseName + "_re");
    auto it_imag = parameterValues.find(baseName + "_im");

    if (it_real != parameterValues.end() && it_imag != parameterValues.end()) {
        // Construct a complex number from the real and imaginary parts
        return std::complex<double>(it_real->second, it_imag->second);
	cout << "get par called correctly for par " << complex<double>(it_real->second, it_imag->second) << endl;
    } else {
        throw std::runtime_error("Error: Real or imaginary part for parameter " + baseName + " not found.");
    }
}



// Setter function: sets the value for a given parameter in the map
void BqDqDqbar::SetParameterValue(const std::string& paramName, double value) {
    parameterValues[paramName] = value;  // Insert or update the value for the given parameter name
}



 //compute decay amplitude for each channel
 void BqDqDqbar::compute_decay_amplitudes(const std::string& channel, const CKMParameters& ckm, bool conjugate) {
   //cout << "computing decay amplitude for channel " << channel << endl;
   amplitudes[channel] = std::complex<double>(0.0, 0.0);

    Parameter amp;
   // Get the CKM elements for the current channel, apply conjugation based on the 'conjugate' flag
    std::complex<double> Vcd = conjugate ? std::conj(ckm.getVcd()) : ckm.getVcd();
    std::complex<double> Vcs = conjugate ? std::conj(ckm.getVcs()) : ckm.getVcs();
    std::complex<double> Vbc = conjugate ? std::conj(ckm.getVcb()) : ckm.getVcb();
    std::complex<double> Vtd = conjugate ? std::conj(ckm.getVtd()) : ckm.getVtd();
    std::complex<double> Vbt = conjugate ? std::conj(ckm.getVtb()) : ckm.getVtb();
    std::complex<double> Vud = conjugate ? std::conj(ckm.getVud()) : ckm.getVud();
    std::complex<double> Vbu = conjugate ? std::conj(ckm.getVub()) : ckm.getVub();
    std::complex<double> Vts = conjugate ? std::conj(ckm.getVts()) : ckm.getVts();
    std::complex<double> Vus = conjugate ? std::conj(ckm.getVus()) : ckm.getVus();
    if (channel == "Bddpdm") {
      amp = Vcd * Vbc * (getPar("E1_dcc_Bddpdm") + getPar("A2_cdc_Bddmdp"))
	- Vud*Vbu * (getPar("P1_GIM_dc_Bddpdm") + getPar("P3_GIM_dc_Bddmdp"))
	- Vtd*Vbt * (getPar("P1_dc_Bddpdm") + getPar("P3_dc_Bddmdp"));
        amplitudes[channel] = amp;
    } else if (channel == "Bddpsdm") {
      amp = Vcs*Vbc*getPar("E1_scc_Bddpsdm") - Vts*Vbt*getPar("P1_sc_Bddpsdm")
	- Vus*Vbu*getPar("P1_GIM_sc_Bddpsdm");
        amplitudes[channel] = amp;
    } else if (channel == "Bpdpd0b") {
      amp = Vcd*Vbc*getPar("E1_dcc_Bpdpd0b") - Vtd*Vbt*getPar("P1_dc_Bpdpd0b")
	+ Vud*Vbu*(getPar("A1_dcu_Bpdpd0b") - getPar("P1_GIM_dc_Bpdpd0b"));
        amplitudes[channel] = amp;
    } else if (channel == "Bpdpsd0b") {
        amp = Vcs*Vbc*getPar("E1_scc_Bpdpsd0b") - Vts*Vbt*getPar("P1_sc_Bpdpsd0b")
	+ Vus*Vbu*(getPar("A1_scu_Bpdpsd0b") - getPar("P1_GIM_sc_Bpdpsd0b"));
        amplitudes[channel] = amp;
    } else if (channel == "Bsdpsdms") {
        amp = Vcs * Vbc * (getPar("E1_scc_Bsdpsdms") + getPar("A2_csc_Bsdmsdps"))
	- Vus*Vbu * (getPar("P1_GIM_sc_Bsdpsdms") + getPar("P3_GIM_sc_Bsdpsdms"))
	- Vts*Vbt * (getPar("P1_sc_Bsdpsdms") + getPar("P3_sc_Bsdpsdms"));
        amplitudes[channel] = amp;
    } else if (channel == "Bsdpdms") {
      amp = Vcd*Vbc*getPar("E1_dcc_Bsdpdms") - Vtd*Vbt*getPar("P1_dc_Bsdpdms")
	- Vud*Vbu*getPar("P1_GIM_dc_Bsdpdms");
        amplitudes[channel] = amp;
    }

 }

   
 //getter for amplitudes
 Parameter BqDqDqbar::get_amplitude(const std::string& channel, const CKMParameters& ckm) {
   compute_decay_amplitudes(channel, ckm, false);  // false for normal
    return amplitudes[channel];
}

 
 //getter for conjugated amplitudes
 Parameter BqDqDqbar::get_conjugate_amplitude(const std::string& channel, const CKMParameters& ckm) {
   compute_decay_amplitudes(channel, ckm, true);  // true for conjugate
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
    double decay_width = (std::pow(G_F, 2) / 2) * (p / (std::pow(m_B, 2) * 16 * pi * h_t)) * (std::norm(amplitude));
   

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
double BqDqDqbar::CalculateC(const Parameter& amplitude, const Parameter& conjugate_amplitude, const std::string& channel, const CKMParameters& ckm) {
    // Get q/p based on whether the channel is Bd or Bs (use parseChannel if necessary)
    auto parsed = parseChannel(channel);
    std::string bMeson = parsed.first;

    std::complex<double> q_p = (bMeson == "Bd") ? ckm.get_q_p_Bd() : ckm.get_q_p_Bs();
    
    // Calculate the ratio lambda = eta * (q/p) * (A_cp / A_conj)
    double eta = cpEigenvalue[channel];  // Assume cpEigenvalue map has eta for each channel
    std::complex<double> lambda = eta * q_p * (amplitude / conjugate_amplitude);

    // Calculate C observable: C = (1 - |lambda|^2) / (1 + |lambda|^2)
    double mod_lambda_squared = std::norm(lambda);
    double C = (1.0 - mod_lambda_squared) / (1.0 + mod_lambda_squared);

    return C;
}

 double BqDqDqbar::CalculateS(const Parameter& amplitude, const Parameter& conjugate_amplitude, const std::string& channel, const CKMParameters& ckm) {
    // Get q/p based on whether the channel is Bd or Bs
    auto parsed = parseChannel(channel);
    std::string bMeson = parsed.first;

    std::complex<double> q_p = (bMeson == "Bd") ? ckm.get_q_p_Bd() : ckm.get_q_p_Bs();
    
    // Calculate the ratio lambda = eta * (q/p) * (A_cp / A_conj)
    double eta = cpEigenvalue[channel];  // Assume cpEigenvalue map has eta for each channel
    std::complex<double> lambda = eta * q_p * (amplitude / conjugate_amplitude);

    // Calculate S observable: S = 2 Im(lambda) / (1 + |lambda|^2)
    double mod_lambda_squared = std::norm(lambda);
    double S = 2.0 * std::imag(lambda) / (1.0 + mod_lambda_squared);

    return S;
}



 std::pair<double, double> BqDqDqbar::CalculatePhiAndLambda(const Parameter& amplitude, const Parameter& conjugate_amplitude, const std::string& channel, const CKMParameters& ckm) {
    // Get q/p based on the B meson (Bs in this case)
    std::complex<double> q_p = ckm.get_q_p_Bs();

    // Compute lambda = (q/p) * (A_cp / A_conj)
    std::complex<double> lambda = q_p * (amplitude / conjugate_amplitude);

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


 //----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
bool isFirstRun = true;  // Declare a flag to check if it's the first run

double BqDqDqbar::LogLikelihood(const std::vector<double>& parameters) {
    std::cout << "LogLikelihood called: " << std::endl;
    int expectedSize = 0;
    for (const auto& channel : channelNames) {
      auto it = channelParameters.find(channel);
      if (it != channelParameters.end()) {
	expectedSize += it->second.size() ;  // Each complex parameter has 2 values (real, imaginary)
      }
    }
    std::cout << "Expected parameters size: " << expectedSize << std::endl;
    if (parameters.size() != expectedSize) {
      std::cerr << "Error: parameters.size() = " << parameters.size() 
              << ", but expectedSize = " << expectedSize << std::endl;
      return -1e30;
    }

    // Ensure that the parameters are valid before using them
    if (parameters.empty()) {
        std::cerr << "Error: Empty parameters vector!" << std::endl;
        return -1e30;
    }

    double ll = 0.0;
   
    // Write to file only on the first run
    if (isFirstRun) {
      //BqDqDqbar::DeclareParameters();  // Call the function to declare parameters
        std::ofstream paramFile("parameterBqDD_indices.txt");
        int index = 0;

        for (const auto& channel : channelNames) {
            auto it = channelParameters.find(channel);
            if (it != channelParameters.end()) {
                for (const auto& paramName : it->second) {
                    paramFile << "Index: " << index << ", Channel: " << channel << ", Parameter: " << paramName << std::endl;
                    index += 2; // For both real and imaginary parts
                }
            } else {
                std::cerr << "Channel " << channel << " not found in channelParameters." << std::endl;
            }
        }
        paramFile.close();
        isFirstRun = false;
    }

    // Unpack the parameters and store them in the parametersValue map
   int index = 0; // This index tracks our position in the 'parameters' vector
   for (const auto& channel : channelNames) { // Loop through each channel name
     auto it = channelParameters.find(channel); // Find parameters for the current channel
     if (it != channelParameters.end()) { // Ensure the channel is found in channelParameters
        // Loop over parameters of this channel, in steps of 2 (since params are stored as real, imaginary pairs)
       for (size_t j = 0; j < it->second.size(); j += 2) {
            // Ensure that we don't go out of bounds when accessing the parameters vector
	 if (index + 1 >= parameters.size()) {
	   std::cerr << "Index out of bounds!" << std::endl;
	   return -1e30;  // Return an error if bounds are exceeded
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

            // Increment the index by 2, because we're handling two values (real, imaginary) at a time
         index += 2;
       }
     } else {
        // If the channel is not found in 'channelParameters', print an error message and return a failure code
       std::cerr << "Channel " << channel << " not found in channelParameters." << std::endl;
       return -1e30; // Return an error code
     }

    

   }

   std::cout << "Finished unpacking parameters. Final index value: " << index << std::endl;


    // Loop over channelNames
    for (const std::string& channel : channelNames) {
      
        // Get the amplitude for the current channel
        Parameter amplitude = get_amplitude(channel, ckm);

        // Get the conjugate amplitude
        Parameter conjugate_amplitude = get_conjugate_amplitude(channel, ckm);
        cout << "finished computing amplitude for channel " << channel << endl;

        //if (isValid(amplitude)) {
            // Compute the predicted branching ratio for this channel
            double br_A = CalculateBR(amplitude, channel);
            double br_conj_A = CalculateBR(conjugate_amplitude, channel);
            double br_predicted = 0.5 * (br_A + br_conj_A);

            // Construct the key for accessing the branching ratio in the meas map
            std::string br_key = "BR" + channel;

            // Get the experimental data for the branching ratio from the meas map
           auto br_it = meas.find(br_key);
           if (br_it != meas.end()) {
	     double br_observed = br_it->second.getMean();
             double br_uncertainty = br_it->second.getSigma();

           // Compute the log-likelihood contribution for this channel (Gaussian likelihood)
             double diff = br_predicted - br_observed;
             ll += -0.5 * (diff * diff) / (br_uncertainty * br_uncertainty);
	     std::cout << "Likelihood contribution for BR_" << channel << " is " << ll << std::endl;

	   } else {
	     std::cerr << "Error: Branching ratio for " << br_key << " not found in meas map." << std::endl;
	   } 

         // Calculate A_CP, C, and S observables if they exist for the channel
           std::string acp_key = "ACP" + channel;
           auto acp_it = meas.find(acp_key);
           if (acp_it != meas.end()) {
	     double acp_predicted = CalculateAcp(amplitude, conjugate_amplitude);
             double acp_observed = acp_it->second.getMean();
             double acp_uncertainty = acp_it->second.getSigma();
	     cout << "the ACP you calculated is " << acp_predicted << " while the experimental value is " << acp_observed << endl;

             double diff_acp = acp_predicted - acp_observed;
             ll += -0.5 * (diff_acp * diff_acp) / (acp_uncertainty * acp_uncertainty);
	     std::cout << "Likelihood contribution for ACP_" << channel << " is " << ll << std::endl;

	   }

            // Continue calculating C and S, as well as log-likelihood
           std::string c_key = "C" + channel;
           std::string s_key = "S" + channel;

           auto c_it = meas.find(c_key);
	   if (c_it != meas.end()) {
	     double c_predicted = CalculateC(amplitude, conjugate_amplitude, channel, ckm);
             double c_observed = c_it->second.getMean();
             double c_uncertainty = c_it->second.getSigma();

	     cout << "the C you calculated is " << c_predicted << " while the experimental value is " << c_observed << endl;

             double diff_c = c_predicted - c_observed;
             ll += -0.5 * (diff_c * diff_c) / (c_uncertainty * c_uncertainty);
	     std::cout << "Likelihood contribution for C_" << channel << " is " << ll << std::endl;

	   }

           auto s_it = meas.find(s_key);
	   if (s_it != meas.end()) {
	     double s_predicted = CalculateS(amplitude, conjugate_amplitude, channel, ckm);
	     double s_observed = s_it->second.getMean();
	     double s_uncertainty = s_it->second.getSigma();

	     cout << "the S you calculated is " << s_predicted << " while the experimental value is " << s_observed << endl;

	     double diff_s = s_predicted - s_observed;
	     ll += -0.5 * (diff_s * diff_s) / (s_uncertainty * s_uncertainty);
	     std::cout << "Likelihood contribution for S_" << channel << " is " << ll << std::endl;

	     
	     //} else {
	     // std::cerr << "Warning: Skipping channel " << channel << " due to missing amplitude." << std::endl;
	     // }
	     } 
    }
    cout << "finished calculating ll from meas, starting corrmeas now..." << endl; 
    cout << "your ll so far from uncorrelated meas is ll = " << ll << endl;
    
    // --------------------------------------------------------
    // Correlated measurements handling using corrmeas map
    for (auto& pair : corrmeas) {
        const auto& corr_key = pair.first;
        auto& corrObservable = pair.second;

        // Extract channels and experiment from the key
        auto channel_and_experiment = extractChannelFromCorrKey(corr_key);
        auto& channels = channel_and_experiment.first;
        auto& experiment = channel_and_experiment.second;
	// Create a prediction vector with the correct size
        int nrows = corrObservable.getObs().GetNrows();
	TVectorD predicted_values(nrows);  // Declare predicted_values here!


       // Ensure that the number of predicted values matches the size of the correlation matrix
       if (nrows != corrObservable.getObs().GetNrows()) {
	 std::cerr << "Error: Size mismatch between predicted values (" << predicted_values.GetNrows() 
                   << ") and correlation matrix rows (" << nrows << ") for key: " << corr_key << std::endl;
	 continue;  // Skip this correlated observable if sizes don't match
       }

       // DEBUG: Print correlation matrix size and predicted_values size
       std::cout << "Correlation matrix rows: " << nrows << std::endl;
       std::cout << "Predicted values size: " << predicted_values.GetNrows() << std::endl;

        // Determine if the key corresponds to ACP or CS observables
        if (corr_key.find("ACP") != std::string::npos) {
            // Loop over channels and calculate ACP for each channel
            for (size_t i = 0; i < channels.size(); ++i) {
                predicted_values[i] = CalculateAcp(get_amplitude(channels[i], ckm), get_conjugate_amplitude(channels[i], ckm));
		//cout << "calculated predicted values for ACP for channel " << channels[i] << " and is " << predicted_values[i] << endl;
            }
        } else if (corr_key.find("CS") != std::string::npos) {
            // Calculate C and S observables for the channel (assumes a single channel for CS observables)
            predicted_values[0] = CalculateC(get_amplitude(channels[0], ckm), get_conjugate_amplitude(channels[0], ckm), channels[0], ckm);
            predicted_values[1] = CalculateS(get_amplitude(channels[0], ckm), get_conjugate_amplitude(channels[0], ckm), channels[0], ckm);
        }/* else if (corr_key.find("phi_lamda") != std::string::npos) {
            // Calculate phi_s and |lambda| observables for the Bsdpsdms channel
            auto phi_lambda_result = CalculatePhiAndLambda(get_amplitude(channels[0], ckm), get_conjugate_amplitude(channels[0], ckm), channels[0], ckm);
            predicted_values[0] = phi_lambda_result.first;  // First index for phi_s
            predicted_values[1] = phi_lambda_result.second; // Second index for |lambda|
	    } */else {
            std::cerr << "Unknown correlation key format: " << corr_key << std::endl;
            continue;
        }

	// Before calling logweight
        if (predicted_values.GetNrows() != corrObservable.getObs().GetNrows()) {
	  std::cerr << "Mismatch in sizes for " << corr_key 
                    << ": predicted_values size " << predicted_values.GetNrows() 
                    << ", correlation matrix size " << corrObservable.getObs().GetNrows() << std::endl;
	  continue;
	}


        // Calculate log-likelihood contribution using correlated data
	cout << "the contribution to ll from correlated observables is " << corrObservable.logweight(predicted_values) << endl;
        ll += corrObservable.logweight(predicted_values);
    }
    // std::cout << "LogLikelihood value: " << ll << std::endl;
    cout << "finished computing corrmeas..." << endl; 
    cout << "your total likelihood is = " << ll << endl; 

    return ll;
}



//---------------------------------------------------
void BqDqDqbar::MCMCUserIterationInterface() {
    std::vector<double> pars;
    
    for (unsigned int i = 0; i < fMCMCNChains; ++i) {
        pars = fMCMCStates.at(i).parameters;
        LogLikelihood(pars);  // Evaluate the likelihood for this parameter set

        int index = 0;  // To loop through the parameters in the correct order

        // Loop over channels and their parameters
        for (const auto& channel : channelNames) {
            auto it = channelParameters.find(channel);
            if (it != channelParameters.end()) {
                for (size_t j = 0; j < it->second.size(); j += 2) {
                    // Real and imaginary parts of the parameter
                    double realPart = pars[index];
                    double imagPart = pars[index + 1];

                    // Fill histograms for real and imaginary parts
		    /*   histos.h1d[channel + "_re"]->Fill(realPart);  // Filling histogram for the real part
                    histos.h1d[channel + "_im"]->Fill(imagPart);  // Filling histogram for the imaginary part

                    // Calculate norm and phase
                    double norm = std::sqrt(realPart * realPart + imagPart * imagPart);
                    double phase = std::atan2(imagPart, realPart);  // Correct phase computation

                    // Fill histograms for norm and phase
                    histos.h1d[channel + "_norm"]->Fill(norm);
                    histos.h1d[channel + "_phase"]->Fill(phase); */

                    index += 2;  // Move to the next pair of parameters
                }
            }
        }
    }
}


//void BqDqDqbar::PrintHistogram() {
//   histos.write();  // Write all histograms to a ROOT file
//}


// ---------------------------------------------------------

    

