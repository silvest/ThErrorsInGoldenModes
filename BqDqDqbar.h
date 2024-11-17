#ifndef BQDQDQBAR_H
#define BQDQDQBAR_H

#include <BAT/BCModel.h>
#include "CKM.h"
#include <BAT/BCH2D.h>
#include <TH2D.h>
#include <map>
#include <random>  
#include <unordered_map>
#include <complex>
#include <string>
#include <tuple>

#include "histo.h"
#include "dato.h"
#include "CorrelatedGaussianObservables.h"

using namespace std;
using Parameter = complex<double>;


class BqDqDqbar : public BCModel {
public:
    // Constructors and destructor
    BqDqDqbar();
    ~BqDqDqbar();
   //map to store all the parameters used in amplitudes
  std::map<std::string, std::vector<std::string>> channelParameters;


    // Declare other member variables
  std::map<std::string, double> parameterValues;  // Stores parameter values
  
 
  
  // map<string,double> obs;
  std::map<std::string, double> obs;



  // Map to store computed amplitudes
    map<string, Parameter> amplitudes;
  
  double LogLikelihood(const std::vector<double> &parameters);
  // void MCMCUserIterationInterface(); 
  void PrintHistogram();

  //vector with all the channel names:
  std::vector<std::string> channelNames = {"Bddpdm", "Bddpsdm", "Bpdpd0b", "Bpdpsd0b", "Bsdpsdms", "Bsdpdms"};


  // Map to store the CP eigenvalue (eta) for each channel
  std::map<std::string, int> cpEigenvalue = {
      {"Bddpdm", +1},  
      {"Bddpsdm", +1}, 
      {"Bsdpsdms", +1},
      {"Bsdpdms", +1}
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
    {"dp", 1.86966}, 
    {"dm", 1.86966}, 
    {"d0b", 1.86484}, 
    {"dps", 1.96835}, 
    {"dms", 1.96835} 
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
  double CalculateC(const Parameter& amplitude, const Parameter& conjugate_amplitude, const std::string& channel);
  double CalculateS(const Parameter& amplitude, const Parameter& conjugate_amplitude, const std::string& channel);
  std::pair<double, double> CalculatePhiAndLambda(const Parameter& amplitude, const Parameter& conjugate_amplitude, const std::string& channel);

  

  std::pair<std::vector<std::string>, std::string> extractChannelFromCorrKey(const std::string& corr_key);

   double Calculate_UncorrelatedObservables(const std::map<std::string, std::pair<Parameter, Parameter>>& amplitude_map);
   double Calculate_CorrelatedObservables(const std::map<std::string, std::pair<Parameter, Parameter>>& amplitude_map);

  void MCMCUserIterationInterface();
  void SaveHistograms(const std::string& filename);
  void PrintObservablePulls(const std::string& filename);



  
private:
    // Existing members
  std::map<std::string, dato> meas;
  std::map<std::string, dato> newmeas;
  histo histos;
  map<string,CorrelatedGaussianObservables> corrmeas;

  
};

#endif
