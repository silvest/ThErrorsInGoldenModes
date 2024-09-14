#include "BqDqDqbar.h"
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

BqDqDqbar::BqDqDqbar() : BCModel(), histos(obs)
{
    //  vector<dato> CorrData;
    //  TMatrixDSym Corr, Corr2;

    d2r = M_PI/180.;
    r2d = (rad ? 1. : 180./M_PI);
  
    std::vector<dato> data;
    PDGAverage pdgaverage;
    double limit1 = 0.0;
    double limit2 = 1000; //all amplitudes calculated from exp BR are of order 1e-5 / G_F, this value is multiplied by 100 accounting for ckm factors


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
    meas.insert(pair<string, dato>("Bpdpd0bVSBpdpsd0b",dato(7.25e-2, 0.09e-2, 0.09e-2 ))); //LHCb:2023wbb


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
    data.push_back(dato(0.265, 0.175, 0.020));  //LHCb:2016inx
    data.push_back(dato(-0.43, 0.16, 0.05)); //Belle:2012mef
    data.push_back(dato(-0.07, 0.23, 0.03));  //BaBar:2008xnt
    data.push_back(dato(-0.91, 0.23, 0.06)); //Belle:2007ebz
    //NUOVA MISURA
    data.push_back(dato(0.162, 0.088, 0.009));  //LHCb:2024gkk
    

    pdgaverage.setData(data);
    pdgaverage.setName("CBddpdm");
    pdgaverage.CalculateAverage();
    
    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    //SBddpdm
    data.push_back(dato(-0.535, 0.165, 0.05));  //LHCb:2016inx
    data.push_back(dato(-0.99, 0.175, 0.08)); //Belle:2012mef
    data.push_back(dato(-0.63, 0.36, 0.05));  //BaBar:2008xnt
    data.push_back(dato(-1.13, 0.37, 0.09)); //Belle:2007ebz
    //NUOVA MISURA
    data.push_back(dato(-0.549, 0.085, 0.015));  //LHCb:2024gkk
    

    pdgaverage.setData(data);
    pdgaverage.setName("SBddpdm");
    pdgaverage.CalculateAverage();
    
    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    //ACPBpdpd0b
    data.push_back(dato(2.5, 1.0, 0.5));  //LHCb:2023wbb
    data.push_back(dato(0, 8, 2)); //Belle:2008doh
    data.push_back(dato(-13, 14, 2));  //BaBar:2006uih
    

    pdgaverage.setData(data);
    pdgaverage.setName("ACPBpdpd0b");
    pdgaverage.CalculateAverage();
    
    meas.insert(pair<string, dato>(pdgaverage.getName(), dato(pdgaverage.getAverage(), pdgaverage.getUncertainty())));
    data.clear();

    //ACPBpdpsd0b
    meas.insert(pair<string, dato>("ACPBpdpsd0b",dato(0.5, 0.2, 0.6))); //LHCb:2023wbb

  //method to define the parameters needed to calculate each decay amplitude, uses BCModel method AddParameter
void BqDqDqbar::DefineParameters(const string& channel) {
  if (channel == "Bddpdm") {
        std::vector<std::string> params = {
            "E1_dcc_Bddpdm",
            "A2_cdc_Bddmdp",
            "P1_GIM_dc_Bddpdm",
            "P3_GIM_dc_Bddmdp",
            "P1_dc_Bddpdm",
            "P3_dc_Bddmdp"
        };
        channelParameters[channel] = params;
        for (const auto& param : params) {
            AddParameter(param, limit1, limit2);
        }
    } else if (channel == "Bddpsdm") {
        std::vector<std::string> params = {
            "E1_scc_Bddpsdm",
            "P1_GIM_sc_Bddpsdm",
            "P1_sc_Bddpsdm",
        channelParameters[channel] = params;
        for (const auto& param : params) {
            AddParameter(param, limit1, limit2);
        }
    } else if (channel == "Bpdpd0b") {
        std::vector<std::string> params = {
            "E1_dcc_Bpdpd0b",
            "P1_dc_Bpdpd0b",
            "A1_dcu_Bpdpd0b",
            "P1_GIM_dc_Bpdpd0b"
        };
        channelParameters[channel] = params;
        for (const auto& param : params) {
            AddParameter(param, limit1, limit2);
        }
    } else if (channel == "Bpdpsd0b") {
        std::vector<std::string> params = {
            "E1_scc_Bpdpsd0b",
            "P1_sc_Bpdpsd0b",
            "A1_scu_Bpdpsd0b",
            "P1_GIM_sc_Bpdpsd0b"
        };
        channelParameters[channel] = params;
        for (const auto& param : params) {
            AddParameter(param, limit1, limit2);
        }
    } else if (channel == "Bsdpsdms") {
        std::vector<std::string> params = {
            "E1_scc_Bsdpsdms",
            "A2_csc_Bsdmsdps",
            "P1_sc_Bsdpsdms",
            "P3_cs_Bsdpsdms",
            "P1_GIM_sc_Bsdpsdms",
            "P3_GIM_sc_Bsdpsdms"
        };
        channelParameters[channel] = params;
        for (const auto& param : params) {
            AddParameter(param, limit1, limit2);
        }
    } else if (channel == "Bsdpdms") {
    std::vector<std::string> params = {
            "E1_dcc_Bsdpsdms",
            "P1_dc_Bsdpsdms",
            "P1_GIM_dc_Bsdpsdms",
    };
    channelParameters[channel] = params;
    for (const auto& param : params) {
        AddParameter(param, limit1, limit2);
    }
  } else {
        // Handle the case where the channel is not recognized
        std::cerr << "Error: Unrecognized channel \"" << channel << "\" in DefineParameters." << std::endl;
        throw std::runtime_error("Unrecognized channel: " + channel);
    }
}

 //function that given the channels and the string inside the map channelParameters adds every parameter name to a map <string, double> and returns said map. now parameterValues[channelParameters] is a complex double. 
std::map<std::string, Parameter> BqDqDqbar::DeclareParameters() {
    std::map<std::string, Parameter> parameterValues;

    // Ensure channelParameters is populated
    for (const auto& channel : channelNames) {
        // Call DefineParameters once per channel
        if (channelParameters.find(channel) == channelParameters.end()) {
            BModel::DefineParameters(channel);
        }

        // Add parameters to parameterValues
        auto it = channelParameters.find(channel);
        if (it != channelParameters.end()) {
            for (const auto& param : it->second) {
	      parameterValues[param] = std::complex<double>(0.0, 0.0); //initialize 
            }
        } else {
            std::cerr << "Channel " << channel << " not found in channelParameters." << std::endl;
        }
    }

    return parameterValues;
}

 // Getter function implementation
Parameter BqDqDqbar::getPar(const std::string& name) const {
    auto it = parameterValues.find(name);
    if (it != parameterValues.end()) {
        return it->second;
    } else {
        throw std::runtime_error("Error: Parameter " + name + " not found.");
    }
}


// Setter function: sets the value for a given parameter in the map
void BqDqDqbar::SetParameterValue(const std::string& paramName, Parameter value) {
    parameterValues[paramName] = value;  // Insert or update the value for the given parameter name
}



 //compute decay amplitude for each channel
 void BqDqDqbar::compute_decay_amplitudes(const std::string& channel, const CKMParameters& ckm, bool conjugate) {
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
BqDqDqbar::~BqDqDqbar()
{
  // default destructor
};

// --------------------------------------------------------- 
 std::pair<std::string, std::pair<std::string, std::string>> BqDqDqbar::parseChannel(const std::string& channel) {
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
    
    // Find two final-state mesons by iterating through the remaining string
    for (const auto& meson : mesonMasses) {
        if (remaining.find(meson.first) == 0) {
            meson1 = meson.first;
            remaining = remaining.substr(meson.first.length());  // Remove meson1 from remaining
            break;
        }
    }

    for (const auto& meson : mesonMasses) {
        if (remaining.find(meson.first) == 0) {
            meson2 = meson.first;
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
    double m1 = mesonMasses[meson1];
    double m2 = mesonMasses[meson2];

    // Calculate the momentum magnitude |p| of the final-state mesons
    double p = sqrt((m_B * m_B - (m1 + m2) * (m1 + m2)) * (m_B * m_B - (m1 - m2) * (m1 - m2))) / (2.0 * m_B);

    // The decay width (Γ) is proportional to |p| and the square of the amplitude
    double decay_width = ((G_F^2) / 2) * (p / (m_B^2 * 16 * pi * h_t)) * (std::norm(amplitude));  // Assuming amplitude is a complex number
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

    // Compute A_CP = (|A|^2 - |Abar|^2) / (|A|^2 + |Abar|^2)
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


 //----------------------------------------------------------------------------------
 bool isFirstRun = true;  // Declare a flag to check if it's the first run

double BqDqDqbar::LogLikelihood(const std::vector<double>& parameters) {
    double ll = 0.0;
    ckm.sampleCKMParameters();  // Re-sample CKM elements at each iteration


    // Write to file only on the first run
    if (isFirstRun) {
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
    int index = 0;
    for (const auto& channel : channelNames) { 
        auto it = channelParameters.find(channel);
        if (it != channelParameters.end()) {
            for (const auto& paramName : it->second) {
                double realPart = parameters[index];
                double imagPart = parameters[index + 1];
                parametersValue[paramName] = std::complex<double>(realPart, imagPart);
                index += 2;
            }
        }
    }
    // Loop over channelNames
    for (const std::string& channel : channelNames) { 
        
        // Get the amplitude for the current channel
      Parameter amplitude = get_amplitude(channel, ckm); 

        // Get the conjugate amplitude
      Parameter conjugate_amplitude = get_conjugate_amplitude(channel, ckm); 

        if (isValid(amplitude)) {

            // Compute the predicted branching ratio for this channel
            double br_A = CalculateBR(amplitude, channel);
            double br_conj_A = CalculateBR(conjugate_amplitude, channel);
            double br_predicted = 0.5 * (br_A + br_conj_A);

            // Construct the key for accessing the branching ratio in the meas map
            std::string br_key = "BR" + channel; 

            // Get the experimental data for the branching ratio from the meas map
            if (meas.find(br_key) != meas.end()) {
                double br_observed = meas[br_key].getMean();
                double br_uncertainty = meas[br_key].getSigma();

                // Compute the log-likelihood contribution for this channel (Gaussian likelihood)
                double diff = br_predicted - br_observed;
                ll += -0.5 * (diff * diff) / (br_uncertainty * br_uncertainty);
            } else {
                std::cerr << "Error: Branching ratio for " << br_key << " not found in meas map." << std::endl;
            }

            // Calculate A_CP, C, and S observables if they exist for the channel
            std::string acp_key = "ACP" + channel; 
            if (meas.find(acp_key) != meas.end()) {
                double acp_predicted = CalculateAcp(amplitude, conjugate_amplitude);
                double acp_observed = meas[acp_key].getMean();
                double acp_uncertainty = meas[acp_key].getSigma();

                double diff_acp = acp_predicted - acp_observed;
                ll += -0.5 * (diff_acp * diff_acp) / (acp_uncertainty * acp_uncertainty);
            }
	   
         // Continue calculating C and S, as well as log-likelihood
            std::string c_key = "C" + channel;
            std::string s_key = "S" + channel;

            if (meas.find(c_key) != meas.end()) {
	      double c_predicted = CalculateC(amplitude, conjugate_amplitude, channel, ckm);
	      double c_observed = meas[c_key].getMean();
	      double c_uncertainty = meas[c_key].getSigma();

	      double diff_c = c_predicted - c_observed;
	      ll += -0.5 * (diff_c * diff_c) / (c_uncertainty * c_uncertainty);
	    }

	    if (meas.find(s_key) != meas.end()) {
	      double s_predicted = CalculateS(amplitude, conjugate_amplitude, channel, ckm);
              double s_observed = meas[s_key].getMean();
              double s_uncertainty = meas[s_key].getSigma();

              double diff_s = s_predicted - s_observed;
              ll += -0.5 * (diff_s * diff_s) / (s_uncertainty * s_uncertainty);
	    }
	} else {
	  std::cerr << "Warning: Skipping channel " << channel << " due to missing amplitude." << std::endl;
	}
    }
    return ll;  // Return total log-likelihood
}


void BModel::MCMCUserIterationInterface()
{
  std::vector<double> pars;
  for (unsigned int i = 0; i < fMCMCNChains; ++i) {
//        std::cout << "interface: chain " << i << ", ";
    pars = fMCMCStates.at(i).parameters;
    //    std::cout << pars.size() << std::endl;
    LogLikelihood(pars);
    // histos.fillh1d();
    // histos.fillh2d();
  }
}

void BModel::PrintHistogram()
{
  histos.write();
}

// ---------------------------------------------------------

    

