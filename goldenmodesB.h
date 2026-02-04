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
    map<string, vector<string>> channelParameters;

    // Declare other member variables
    map<string, double> parameterValues; // Stores parameter values

    // map<string,double> obs;
    map<string, double> obs;

    double LogLikelihood(const vector<double> &parameters);
    // void MCMCUserIterationInterface();
    void PrintHistogram();

    // vector with all the channel names:
    vector<string> channelNames = {"Bpjpsipp", "Bpjpsikp", "Bdjpsip0", "Bdjpsik0s", "Bdjpsik0l", "Bdjpsik0", "Bdjpsieta", "Bdjpsietap", "Bsjpsip0", "Bsjpsik0s", "Bsjpsik0l", "Bsjpsik0b", "Bsjpsieta", "Bsjpsietap", "Bpjpsirp", "Bpjpsikstp", "Bdjpsirho0", "Bdjpsikst0", "Bdjpsiphi", "Bdjpsiom", "Bsjpsirho0", "Bsjpsikbst0", "Bsjpsiphi", "Bsjpsiom", "Bpdpd0b", "Bpdspd0b", "Bddpdm", "Bddspdm", "Bddspdsm", "Bdd0d0b", "Bsdpdm", "Bsdpdsm", "Bsdspdsm", "Bsd0d0b"};
    vector<string> channelNamesSU3;
    vector<string> pseudoscalarMesonChannelsSU3 = {"Bdjpsik0", "Bdjpsip0", "Bdjpsieta8", "Bdjpsieta1", "Bpjpsikp", "Bpjpsipp", "Bsjpsip0", "Bsjpsik0b", "Bsjpsieta8", "Bsjpsieta1"};
    vector<string> vectorMesonChannels = {
        "Bpjpsirhop", "Bpjpsikstp", "Bdjpsirho0", "Bdjpsikst0", "Bdjpsiphi", "Bdjpsiom", "Bsjpsirho0", "Bsjpsikbst0", "Bsjpsiphi", "Bsjpsiom"};
    vector<string> vectorMesonChannelsSU3 = {"Bsjpsiphi", "Bsjpsiom", "Bsjpsikbst0", "Bsjpsirho0", "Bdjpsiom", "Bdjpsikst0", "Bdjpsirho0", "Bdjpsiphi", "Bpjpsikstp", "Bpjpsirhop"};
    vector<string> pseudoscalarMesonChannels = {"Bpjpsipp", "Bpjpsikp", "Bdjpsip0", "Bdjpsik0s", "Bdjpsik0l", "Bdjpsik0", "Bdjpsieta", "Bdjpsietap", "Bsjpsip0", "Bsjpsik0s", "Bsjpsik0l", "Bsjpsik0b", "Bsjpsieta", "Bsjpsietap"};
    vector<string> ddbarChannels = {"Bpdpd0b", "Bpdspd0b", "Bddpdm", "Bddspdm", "Bddspdsm", "Bdd0d0b", "Bsdpdm", "Bsdpdsm", "Bsdspdsm", "Bsd0d0b"};
    vector<string> ddbarChannelsSU3 = {"Bsdspdsm", "Bsdpdsm", "Bsdpdm", "Bsd0d0b", "Bddspdsm", "Bddspdm", "Bddpdm", "Bdd0d0b", "Bpdpd0b", "Bpdspd0b"};

    vector<string> channels;

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
    map<string, double> mesonMasses = {
        {"jpsi", 3.09690},
        {"k0", 0.49761},
        {"k0s", 0.49761},
        {"k0l", 0.49761},
        {"k0b", 0.49761}, // K̄⁰ meson
        {"p0", 0.134976},
        {"kp", 0.49367},
        {"pp", 0.139570},
        {"om", 0.78266},
        {"ph", 1.019461},
        {"rh", 1.465},
        {"rho0", 0.77526}, // ρ⁰ meson
        {"rhop", 0.77511}, // ρ⁺ meson
        {"kst", 0.845},
        {"kbst", 0.845},   // K̄* meson
        {"dp", 1.86966},   // D⁺ meson
        {"dm", 1.86966},   // D⁻ meson
        {"d0b", 1.86484},  // D̄⁰ meson
        {"d0", 1.86484},   // D⁰ meson
        {"dsp", 1.96835},  // Ds⁺ meson
        {"dsm", 1.96835},  // Ds⁻ meson
        {"eta", 0.547862}, // eta meson
        {"etap", 0.95778}  // eta' meson
    };

    // to split up the channel name
    pair<string, pair<string, string>> parseChannel(const string &channel) const;

    // function to get tau of B meson
    double getBMesonLifetime(const string &bMeson) const;

    // function to get mass of the B meson
    double getBMesonMass(const string &bMeson) const;

    // define parameters for each channel with AddParameter
    void DefineParameters(const string &channel);

    // declare parameters as double, default initialization
    map<string, double> DeclareParameters();

    // getter for parameterValues (WARNING: you first need to call declareParameters)
    TComplex getPar(const string &name) const;

    // Setter for parameterValues map
    void SetParameterValue(const string &paramName, double value);

    // to get just real or imaginary parts of effective parameters
    double getParameterValue(const string &paramName) const;

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

    map<string, double> CalculatePolarizations(
        pair<TComplex, TComplex> &amplitude_0, pair<TComplex, TComplex> &amplitude_paral, pair<TComplex, TComplex> &amplitude_perp);

    // Function to compute decay amplitudes
    void compute_decay_amplitudes(const string &channel);

    // calculate observables
    double CalculateBR(TComplex amplitude, TComplex amplitude_conj, const string &channel) const;
    double CalculateAcp(const TComplex &amplitude, const TComplex &conjugate_amplitude) const;
    double CalculateAlpha(const TComplex &amplitude, const TComplex &conjugate_amplitude, const string &channel);
    double CalculateC(const TComplex &amplitude, const TComplex &conjugate_amplitude, const string &channel);
    pair<double,double> CalculateS(const TComplex &amplitude, const TComplex &conjugate_amplitude, const string &channel);
    tuple<double, double, double> CalculatePhiAndLambda(const TComplex &amplitude, const TComplex &conjugate_amplitude, const string &channel);
    pair<vector<string>, string> extractChannelFromCorrKey(const string &corr_key);
    map<string, double> getPolarizationParams(const string &channel, const map<string, pair<TComplex, TComplex>> &amplitude_map);

    double Calculate_UncorrelatedObservables(map<string, pair<TComplex, TComplex>> &amplitude_map);
    double Calculate_CorrelatedObservables(map<string, pair<TComplex, TComplex>> &amplitude_map);

    void MCMCUserIterationInterface();
    void SaveHistograms(const string &filename);
    void PrintObservablePulls(const string &filename);

private:
    TComplex lam_bs_c;
    TComplex lam_bs_u;
    TComplex lam_bd_c;
    TComplex lam_bd_u;
    TComplex lamst_bs_c;
    TComplex lamst_bs_u;
    TComplex lamst_bd_c;
    TComplex lamst_bd_u;

    string addPolarizationSuffix(string amplitude, string suffix) const
    {
        size_t pos_re = amplitude.find("_re");
        size_t pos_im = amplitude.find("_im");

        if (pos_re != string::npos)
        {
            amplitude.insert(pos_re, suffix);
        }
        else if (pos_im != string::npos)
        {
            amplitude.insert(pos_im, suffix);
        }
        else
        {
            cerr << "Target substring not found!" << endl;
            exit(1);
        }
        return amplitude;
    }

    int evaluatedevts;

    void addAmplitudeParameter(const string &amplitudeName, const double &lowerLimit, const double &upperLimit, bool isVectorChannel = false)
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
            vector<string> polarizations = {"_0", "_paral", "_perp"};
            for (const auto &pol : polarizations)
            {
                string fullAmplitudeName = addPolarizationSuffix(amplitudeName, pol);
                if (!referenceAmplitudes.count(fullAmplitudeName))
                {
                    AddParameter(fullAmplitudeName, lowerLimit, upperLimit);
                    referenceAmplitudes.emplace(fullAmplitudeName);
                }
                // If already defined, just skip (shared amplitude across channels)
            }
        }
    }

    void addSU3BreakingParameter(const string &amplitudeName, const string &baseParameterName, bool isVectorChannel = false)
    {
        // For vector channels, check for the _0 polarization suffix since that's how they're stored
        string checkName = isVectorChannel ? addPolarizationSuffix(baseParameterName, "_0") : baseParameterName;
        // Remember we store the reference amplitude in deltaReferenceAmplitudes map without the _re/_im suffix
        string checkBaseName = baseParameterName;;
        size_t pos;
        if ((pos = checkBaseName.find("_re")) != string::npos ||
            (pos = checkBaseName.find("_im")) != string::npos)
        {
            checkBaseName = checkBaseName.substr(0, pos);
        }
        checkBaseName = isVectorChannel ? checkBaseName + "_0" : checkBaseName;
        if (referenceAmplitudes.count(checkName) || deltaReferenceAmplitudes.find(checkBaseName) != deltaReferenceAmplitudes.end())
        {
            if (dsu3_limit > 0.0)
            {
                if (amplitudeName.find("EW") != string::npos)
                {
                    // add EWP SU(3) breaking parameters only if EWP contributions are nonvanishing
                    if (ewp_limit > 0.0)
                        addAmplitudeParameter(amplitudeName, -dsu3_limit, dsu3_limit, isVectorChannel);
                    else
                        addAmplitudeParameter(amplitudeName, 0.0, 0.0, isVectorChannel); // fix to zero
                }
                else
                    addAmplitudeParameter(amplitudeName, -dsu3_limit, dsu3_limit, isVectorChannel);

                //Register the breaking amplitude as available for chaining
                string deltaBase = amplitudeName;
                string refBase = baseParameterName;
                size_t pos;
                // Remove "delta_" prefix for mapping
                if ((pos = deltaBase.find("delta_")) != string::npos)
                {
                    deltaBase = deltaBase.substr(pos + 6); // length of "delta_" is 6
                }
                if ((pos = deltaBase.find("_re")) != string::npos ||
                    (pos = deltaBase.find("_im")) != string::npos)
                {
                    deltaBase = deltaBase.substr(0, pos);
                }
                if ((pos = refBase.find("_re")) != string::npos || 
                    (pos = refBase.find("_im")) != string::npos)
                {
                    refBase = refBase.substr(0, pos);
                }

                if (!isVectorChannel)
                {
                    deltaReferenceAmplitudes[deltaBase]= refBase;
                }
                else
                {
                    // For vector channels, register all polarization variants
                    vector<string> polarizations = {"_0", "_paral", "_perp"};
                    for (const auto &pol : polarizations)
                    {
                        deltaReferenceAmplitudes[deltaBase + pol] = refBase + pol;
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
    map<string, dato> meas;
    map<string, dato> newmeas;
    Histos histos;
    map<string, CorrelatedGaussianObservables> corrmeas;
    map<string, vector<string>> corrmeas_channels;
    // vector to store reference amplitudes for SU(3) breaking
    unordered_set<string> referenceAmplitudes;
    map<string, pair<TComplex, TComplex>> amplitude_map;
    // Map from delta parameter to its reference parameter (for SU(3) breaking chain)
    // e.g., "delta_E2t_ccdd_BJPSIP" -> "E2t_ccsd_BJPSIP"
    map<string, string> deltaReferenceAmplitudes;
};

#endif
