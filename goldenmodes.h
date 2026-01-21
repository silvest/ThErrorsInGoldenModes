#ifndef GOLDENMODES_H
#define GOLDENMODES_H

#include <BAT/BCModel.h>
#include "CKM.h"
#include <BAT/BCH2D.h>
#include <TH2D.h>
#include <map>
#include <random>
#include <unordered_map>
#include <TComplex.h>
#include <string>
#include <tuple>
#include <algorithm>


#include "histo.h"
#include "dato.h"
#include "CorrelatedGaussianObservables.h"

using namespace std;
using Parameter = complex<double>;


class goldenmodes : public BCModel {
public:
    // Constructors and destructor
    goldenmodes();
    ~goldenmodes();
   //map to store all the parameters used in amplitudes
  std::map<std::string, std::vector<std::string>> channelParameters;


    // Declare other member variables
  std::map<std::string, double> parameterValues;  // Stores parameter values



  // map<string,double> obs;
  std::map<std::string, double> obs;



  // Map to store computed amplitudes
    map<string, Parameter> amplitudes;
    map<string, Parameter> polarized_amp;

  double LogLikelihood(const std::vector<double> &parameters);
  // void MCMCUserIterationInterface();
  void PrintHistogram();

  //vector with all the channel names:
  std::vector<std::string> channelNames = {"Bpjpsipp", "Bpjpsikp", "Bdjpsip0", "Bdjpsik0s", "Bdjpsik0l", "Bdjpsik0", "Bdjpsiet", "Bdjpsietp", "Bsjpsip0", "Bsjpsik0s", "Bsjpsik0l", "Bsjpsik0b", "Bsjpsiet", "Bsjpsietp", "Bpjpsirp", "Bpjpsikstp", "Bdjpsir0", "Bdjpsikst0", "Bdjpsiph", "Bdjpsiom", "Bsjpsir0", "Bsjpsikst0", "Bsjpsiph", "Bsjpsiom", "Bpdpd0b", "Bpdspd0b", "Bddpdm", "Bddpsdm", "Bddpsdms", "Bdd0d0b", "Bsdpdm", "Bsdpdms", "Bsdpsdms", "Bsd0d0b"};
  std::vector<std::string> channelNamesSU3 = {"Bpjpsipp", "Bpjpsikp", "Bdjpsip0", "Bdjpsik0", "Bdjpset1", "Bdjpsiet8", "Bsjpsip0", "Bsjpsik0b", "Bsjpsiet1", "Bsjpsiet8", "Bpjpsirp", "Bpjpsikstp", "Bdjpsir0", "Bdjpsikst0", "Bdjpsiph", "Bdjpsiom", "Bsjpsir0", "Bsjpsikst0", "Bsjpsiph", "Bsjpsiom", "Bpdpd0b", "Bpdspd0b", "Bddpdm", "Bddpsdm", "Bddpsdms", "Bdd0d0b", "Bsdpdm", "Bsdpdms", "Bsdpsdms", "Bsd0d0b"};
  std::vector<std::string> vectorMesonChannels = {
   "Bpjpsirp", "Bpjpsikstp", "Bdjpsir0", "Bdjpsikst0", "Bdjpsiph", "Bdjpsiom", "Bsjpsir0", "Bsjpsikst0", "Bsjpsiph", "Bsjpsiom"
};



  // Map to store the CP eigenvalue (eta) for each channel
  std::map<std::string, int> cpEigenvalue = {
      {"Bdjpsik0s", +1},
      {"Bdjpsik0l", -1},
      {"Bdjpsip0", +1},
      {"Bdjpsiom", +1},
      {"Bsjpsik0s", +1},
      {"Bdjpsikst", +1},
      {"Bdjpsirh", +1},
      {"Bddpdm", +1},      // Bd→D⁺D⁻: two pseudoscalars, CP-even
      {"Bsdpsdms", +1}     // Bs→Ds⁺Ds⁻: two pseudoscalars, CP-even
};

  //global variables
  const double m_Bp = 5.27941;  // B+ meson mass in GeV
  const double m_Bd = 5.27972;    // B^0_d meson mass in GeV
  const double m_Bs = 5.36693;    // B^0_s meson mass in GeV

  const double tau_Bp = 1.638e-12;  // B+ lifetime in seconds
  const double tau_Bd = 1.517e-12;    // B^0_d lifetime in seconds
  const double tau_Bs = 1.5195e-12;    // B^0_s lifetime in seconds

  const double h_t = 6.582119569e-25; //h tagliato in Gev s
  const double G_F = 1.1663788e-5; //Fermi constant in  GeV



  //map with the final state meson masses:
  std::map<std::string, double> mesonMasses = {
    {"jpsi", 3.09690},
    {"k0", 0.49761},
    {"k0s", 0.49761},
    {"k0l", 0.49761},
    {"p0", 0.134976},
    {"kp", 0.49367},
    {"pp", 0.139570},
    {"om", 0.78266},
    {"ph", 1.019461},
    {"rh", 1.465},
    {"kst", 0.845},
    {"dp", 1.86966},    // D⁺ meson
    {"dm", 1.86966},    // D⁻ meson
    {"d0b", 1.86484},   // D̄⁰ meson
    {"dps", 1.96835},   // Ds⁺ meson
    {"dms", 1.96835}    // Ds⁻ meson
};


  //to split up the channel name
  std::pair<std::string, std::pair<std::string, std::string>> parseChannel(const std::string& channel) const;

  //function to get tau of B meson
  double getBMesonLifetime(const std::string& bMeson) const;

  //function to get mass of the B meson
  double getBMesonMass(const std::string& bMeson) const;

  //define parameters for each channel with AddParameter
  void DefineParameters(const string& channel);

  //declare parameters as double, default initialization
  std::map<std::string, double> DeclareParameters();

  //getter for parameterValues (WARNING: you first need to call declareParameters)
  Parameter getPar(const std::string& name) const;


  // Setter for parameterValues map
  void SetParameterValue(const std::string& paramName, double value);


  //to get just real or imaginary parts of effective parameters
  double getParameterValue(const std::string& paramName) const;


  // getter for meas, the map with exp measures
  map<string,dato> getMeas() const {
      return meas;
  }

  // setter for meas
  void setMeas(const map<string,dato>& newMeas) {
    meas = newMeas;
  }

    // Utility methods

    // Function to compute decay amplitudes
  void compute_decay_amplitudes(const std::string& channel, bool conjugate);

    // Accessor for specific channel amplitude
  Parameter get_amplitude(const std::string& channel);
  Parameter get_conjugate_amplitude(const std::string& channel);


  //calculate observables
  double CalculateBR(Parameter amplitude, const string& channel) const;
  double CalculateAcp(const Parameter& amplitude, const Parameter& conjugate_amplitude) const;
  double CalculateAlpha(const Parameter& amplitude, const Parameter& conjugate_amplitude, const std::string& channel);
  double CalculateC(const Parameter& amplitude, const Parameter& conjugate_amplitude, const std::string& channel);
  double CalculateS(const Parameter& amplitude, const Parameter& conjugate_amplitude, const std::string& channel);
  std::pair<double, double> CalculatePhiAndLambda(const Parameter& amplitude, const Parameter& conjugate_amplitude, const std::string& channel);

  std::pair<std::vector<std::string>, std::string> extractChannelFromCorrKey(const std::string& corr_key);
  map<string, double> getPolarizationParams(const string& channel, const std::map<std::string, std::pair<Parameter, Parameter>>& amplitude_map);

  double Calculate_UncorrelatedObservables(const std::map<std::string, std::pair<Parameter, Parameter>>& amplitude_map);
  double Calculate_CorrelatedObservables(const std::map<std::string, std::pair<Parameter, Parameter>>& amplitude_map);

  void MCMCUserIterationInterface();
  void SaveHistograms(const std::string& filename);
  void PrintObservablePulls(const std::string& filename);




private:
    // Existing members
  std::map<std::string, dato> meas;
  std::map<std::string, dato> newmeas;
  Histos histos;
  map<string,CorrelatedGaussianObservables> corrmeas;


};

#endif
