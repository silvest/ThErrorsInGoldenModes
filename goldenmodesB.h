#ifndef GOLDENMODESB_H
#define GOLDENMODESB_H

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
#include <unordered_set>

#include "histo.h"
#include "dato.h"
#include "CorrelatedGaussianObservables.h"

using namespace std;

class goldenmodesB : public BCModel
{
public:
    // Constructors and destructor
    goldenmodesB(double &dsu3_limit, double &ewp_limit, bool BJPSIP = true, bool BJPSIV = true, bool BDDb = true);
    ~goldenmodesB();
    // map to store all the parameters used in amplitudes
    std::map<std::string, std::vector<std::string>> channelParameters;

    // Declare other member variables
    std::map<std::string, double> parameterValues; // Stores parameter values

    // map<string,double> obs;
    std::map<std::string, double> obs;

    // Map to store computed amplitudes
    map<string, TComplex> amplitudes;
    map<string, TComplex> polarized_amp;

    double LogLikelihood(const std::vector<double> &parameters);
    // void MCMCUserIterationInterface();
    void PrintHistogram();

    // vector with all the channel names:
    std::vector<std::string> channelNames = {"Bpjpsipp", "Bpjpsikp", "Bdjpsip0", "Bdjpsik0s", "Bdjpsik0l", "Bdjpsik0", "Bdjpsieta", "Bdjpsietap", "Bsjpsip0", "Bsjpsik0s", "Bsjpsik0l", "Bsjpsik0b", "Bsjpsieta", "Bsjpsietap", "Bpjpsirp", "Bpjpsikstp", "Bdjpsirho0", "Bdjpsikst0", "Bdjpsiphi", "Bdjpsiom", "Bsjpsirho0", "Bsjpsikbst0", "Bsjpsiphi", "Bsjpsiom", "Bpdpd0b", "Bpdspd0b", "Bddpdm", "Bddspdm", "Bddspdsm", "Bdd0d0b", "Bsdpdm", "Bsdpdsm", "Bsdspdsm", "Bsd0d0b"};
    std::vector<std::string> channelNamesSU3 = {"Bdjpsik0", "Bdjpsip0", "Bdjpsieta8", "Bdjpsieta1", "Bpjpsikp", "Bpjpsipp", "Bsjpsip0", "Bsjpsik0b", "Bsjpsieta8", "Bsjpsieta1", "Bsjpsiphi", "Bsjpsiom", "Bsjpsikbst0", "Bsjpsirho0", "Bdjpsiom", "Bdjpsikst0", "Bdjpsirho0", "Bdjpsiphi", "Bpjpsikstp", "Bpjpsirhop", "Bsdspdsm", "Bsdpdsm", "Bsdpdm", "Bsd0d0b", "Bddspdsm", "Bddspdm", "Bddpdm", "Bdd0d0b", "Bpdpd0b", "Bpdspd0b"};
    std::vector<std::string> vectorMesonChannels = {
        "Bpjpsirhop", "Bpjpsikstp", "Bdjpsirho0", "Bdjpsikst0", "Bdjpsiphi", "Bdjpsiom", "Bsjpsirho0", "Bsjpsikbst0", "Bsjpsiphi", "Bsjpsiom"};
    std::vector<std::string> pseudoscalarMesonChannels = {"Bpjpsipp", "Bpjpsikp", "Bdjpsip0", "Bdjpsik0s", "Bdjpsik0l", "Bdjpsik0", "Bdjpsieta", "Bdjpsietap", "Bsjpsip0", "Bsjpsik0s", "Bsjpsik0l", "Bsjpsik0b", "Bsjpsieta", "Bsjpsietap"};
    std::vector<std::string> ddbarChannels = {"Bpdpd0b", "Bpdspd0b", "Bddpdm", "Bddspdm", "Bddspdsm", "Bdd0d0b", "Bsdpdm", "Bsdpdsm", "Bsdspdsm", "Bsd0d0b"};

    std::vector<std::string> channels;

    // global variables
    const double m_Bp = 5.27941; // B+ meson mass in GeV
    const double m_Bd = 5.27972; // B^0_d meson mass in GeV
    const double m_Bs = 5.36693; // B^0_s meson mass in GeV

    const double tau_Bp = 1.638e-12;  // B+ lifetime in seconds
    const double tau_Bd = 1.517e-12;  // B^0_d lifetime in seconds
    const double tau_Bs = 1.5195e-12; // B^0_s lifetime in seconds

    const double h_t = 6.582119569e-25; // h tagliato in Gev s
    const double G_F = 1.1663788e-5;    // Fermi constant in  GeV

    // map with the final state meson masses:
    std::map<std::string, double> mesonMasses = {
        {"jpsi", 3.09690},
        {"k0", 0.49761},
        {"k0s", 0.49761},
        {"k0l", 0.49761},
        {"k0b", 0.49761},  // K̄⁰ meson
        {"p0", 0.134976},
        {"kp", 0.49367},
        {"pp", 0.139570},
        {"om", 0.78266},
        {"ph", 1.019461},
        {"rh", 1.465},
        {"rho0", 0.77526},  // ρ⁰ meson
        {"rhop", 0.77511},  // ρ⁺ meson
        {"kst", 0.845},
        {"kbst", 0.845},    // K̄* meson
        {"dp", 1.86966},  // D⁺ meson
        {"dm", 1.86966},  // D⁻ meson
        {"d0b", 1.86484}, // D̄⁰ meson
        {"d0", 1.86484},  // D⁰ meson
        {"dsp", 1.96835}, // Ds⁺ meson
        {"dsm", 1.96835},  // Ds⁻ meson
        {"eta", 0.547862}, // eta meson
        {"etap", 0.95778} // eta' meson
    };

    // to split up the channel name
    std::pair<std::string, std::pair<std::string, std::string>> parseChannel(const std::string &channel) const;

    // function to get tau of B meson
    double getBMesonLifetime(const std::string &bMeson) const;

    // function to get mass of the B meson
    double getBMesonMass(const std::string &bMeson) const;

    // define parameters for each channel with AddParameter
    void DefineParameters(const string &channel);

    // declare parameters as double, default initialization
    std::map<std::string, double> DeclareParameters();

    // getter for parameterValues (WARNING: you first need to call declareParameters)
    TComplex getPar(const std::string &name) const;

    // Setter for parameterValues map
    void SetParameterValue(const std::string &paramName, double value);

    // to get just real or imaginary parts of effective parameters
    double getParameterValue(const std::string &paramName) const;

    // getter for meas, the map with exp measures
    map<string, dato> getMeas() const
    {
        return meas;
    }

    // setter for meas
    void setMeas(const map<string, dato> &newMeas)
    {
        meas = newMeas;
    }

    // Utility methods

    std::map<std::string, double> CalculatePolarizations(
std::pair<TComplex, TComplex> &amplitude_0, std::pair<TComplex, TComplex> &amplitude_paral, std::pair<TComplex, TComplex> &amplitude_perp);

    // Function to compute decay amplitudes
    void compute_decay_amplitudes(const std::string &channel, bool conjugate);

    // Accessor for specific channel amplitude
    TComplex get_amplitude(const std::string &channel);
    TComplex get_conjugate_amplitude(const std::string &channel);

    // calculate observables
    double CalculateBR(TComplex amplitude, TComplex amplitude_conj, const string &channel) const;
    double CalculateAcp(const TComplex &amplitude, const TComplex &conjugate_amplitude) const;
    double CalculateAlpha(const TComplex &amplitude, const TComplex &conjugate_amplitude, const std::string &channel);
    double CalculateC(const TComplex &amplitude, const TComplex &conjugate_amplitude, const std::string &channel);
    double CalculateS(const TComplex &amplitude, const TComplex &conjugate_amplitude, const std::string &channel);
    std::pair<double, double> CalculatePhiAndLambda(const TComplex &amplitude, const TComplex &conjugate_amplitude, const std::string &channel);
    std::pair<std::vector<std::string>, std::string> extractChannelFromCorrKey(const std::string &corr_key);
    map<string, double> getPolarizationParams(const string &channel, const std::map<std::string, std::pair<TComplex, TComplex>> &amplitude_map);

    double Calculate_UncorrelatedObservables(std::map<std::string, std::pair<TComplex, TComplex>> &amplitude_map);
    double Calculate_CorrelatedObservables(std::map<std::string, std::pair<TComplex, TComplex>> &amplitude_map);

    void MCMCUserIterationInterface();
    void SaveHistograms(const std::string &filename);
    void PrintObservablePulls(const std::string &filename);

private:
    std::string addPolarizationSuffix(std::string amplitude, std::string suffix) const
    {
        size_t pos_re = amplitude.find("_re");
        size_t pos_im = amplitude.find("_im");

        if (pos_re != std::string::npos)
        {
            amplitude.insert(pos_re, suffix);
        }
        else if (pos_im != std::string::npos)
        {
            amplitude.insert(pos_im, suffix);
        }
        else
        {
            std::cerr << "Target substring not found!" << std::endl;
            exit(1);
        }
        return amplitude;
    }

    void addAmplitudeParameter(const std::string &amplitudeName, const double &lowerLimit, const double &upperLimit, bool isVectorChannel = false)
    {
        if (!isVectorChannel)
        {
            if (!referenceAmplitudes.count(amplitudeName))
            {
                AddParameter(amplitudeName, lowerLimit, upperLimit);
                referenceAmplitudes.emplace(amplitudeName);
            }
            // If already defined, just skip (shared amplitude across channels)
        }
        else
        {
            // For vector channels, define three polarization amplitudes
            std::vector<std::string> polarizations = {"_0", "_paral", "_perp"};
            for (const auto &pol : polarizations)
            {
                std::string fullAmplitudeName = addPolarizationSuffix(amplitudeName, pol);
                if (!referenceAmplitudes.count(fullAmplitudeName))
                {
                    AddParameter(fullAmplitudeName, lowerLimit, upperLimit);
                    referenceAmplitudes.emplace(fullAmplitudeName);
                }
                // If already defined, just skip (shared amplitude across channels)
            }
        }
    }

    void addSU3BreakingParameter(const std::string &amplitudeName, const std::string &baseParameterName, bool isVectorChannel = false)
    {
        // For vector channels, check for the _0 polarization suffix since that's how they're stored
        std::string checkName = isVectorChannel ? addPolarizationSuffix(baseParameterName, "_0") : baseParameterName;
        if (referenceAmplitudes.count(checkName))
        {
            if(amplitudeName.find("EW") != std::string::npos)
            {
                // add EWP SU(3) breaking parameters only if EWP contributions are nonvanishing
                if (ewp_limit > 0.0)
                    addAmplitudeParameter(amplitudeName, -dsu3_limit, dsu3_limit, isVectorChannel);
                else
                    addAmplitudeParameter(amplitudeName, 0.0, 0.0, isVectorChannel); // fix to zero
            }
            else
                addAmplitudeParameter(amplitudeName, -dsu3_limit, dsu3_limit, isVectorChannel);
            
            // Register the derived amplitude (without "delta_" prefix) as available for chaining
            // e.g., if we added delta_E2t_ccsd, then E2t_ccsd is now available as a reference
            if (amplitudeName.rfind("delta_", 0) == 0)
            {
                std::string derivedAmplitude = amplitudeName.substr(6); // Remove "delta_" prefix
                
                // Store the mapping: derived amplitude base name -> reference amplitude base name
                // Extract base names (remove _re/_im suffix)
                std::string deltaBase = amplitudeName;
                std::string refBase = baseParameterName;
                size_t pos;
                if ((pos = deltaBase.find("_re")) != std::string::npos || 
                    (pos = deltaBase.find("_im")) != std::string::npos)
                {
                    deltaBase = deltaBase.substr(0, pos);
                }
                if ((pos = refBase.find("_re")) != std::string::npos || 
                    (pos = refBase.find("_im")) != std::string::npos)
                {
                    refBase = refBase.substr(0, pos);
                }
                // Store: derivedAmplitude (without delta_) -> refBase
                // e.g., "E2t_ccdd_BJPSIP" -> "E2t_ccsd_BJPSIP"
                std::string derivedBase = deltaBase.substr(6); // Remove "delta_" from deltaBase
                if (!isVectorChannel)
                {
                    deltaToReference[derivedBase] = refBase;
                    referenceAmplitudes.emplace(derivedAmplitude);
                }
                else
                {
                    // For vector channels, register all polarization variants
                    std::vector<std::string> polarizations = {"_0", "_paral", "_perp"};
                    for (const auto &pol : polarizations)
                    {
                        // Store mapping for each polarization
                        deltaToReference[derivedBase + pol] = refBase + pol;
                        referenceAmplitudes.emplace(addPolarizationSuffix(derivedAmplitude, pol));
                    }
                }
            }
        }
        else
        {
            cerr << "Error: Reference amplitude " << baseParameterName << " not defined so cannot define SU(3)-related parameters." << endl;
            exit(1);
        }
    }

    double dsu3_limit = 0.0, ewp_limit = 0.0;
    // map to store experimental measurements
    std::map<std::string, dato> meas;
    std::map<std::string, dato> newmeas;
    Histos histos;
    map<string, CorrelatedGaussianObservables> corrmeas;
    map<string, std::vector<std::string>> corrmeas_channels;
    // vector to store reference amplitudes for SU(3) breaking
    std::unordered_set<std::string> referenceAmplitudes;
    std::map<std::string, std::pair<TComplex, TComplex>> amplitude_map;
    // Map from delta parameter to its reference parameter (for SU(3) breaking chain)
    // e.g., "delta_E2t_ccdd_BJPSIP" -> "E2t_ccsd_BJPSIP"
    std::map<std::string, std::string> deltaToReference;
};

#endif
