#ifndef GOLDENMODESB_INDSU3_H
#define GOLDENMODESB_INDSU3_H

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

// -------------------------------------------------------------------
// goldenmodesB_indSU3: Independent-SU3 variant of goldenmodesB
//
// All tilde amplitude parameters (Ẽ₂, P̃₂^GIM, ẼA₂, Ẽ₁, etc.) that
// were previously related by the delta-chain SU(3) breaking scheme are
// now treated as fully independent BAT parameters.
//
// SU(3) symmetry is imposed as a soft Gaussian constraint in the
// log-likelihood for each pair of adjacent SU(3)-related amplitudes:
//   ΔlogL = -0.5 * |A₁ - A₂|² / ( ((|A₁|+|A₂|)/2)² * σ_SU3² )
//
// EW penguin parameters carry a Gaussian prior with σ = ewp_limit:
//   ΔlogL = -0.5 * |A_EWP|² / ewp_limit²
//
// Constructor arguments:
//   ewp_limit   : σ of the Gaussian prior on EW penguin amplitudes
//                 (0 → EW penguins fixed to zero)
//   BJPSIP/V/BDDb : which channel classes to include
//   su3_sigma_in : if > 0, used as fixed σ_SU3
//                  if ≤ 0, σ_SU3 is added as a free BAT parameter
// -------------------------------------------------------------------
class goldenmodesB_indSU3 : public BCModel
{
public:
    goldenmodesB_indSU3(double &ewp_limit, bool BJPSIP = true, bool BJPSIV = true, bool BDDb = true, double su3_sigma_in = 0.3);
    ~goldenmodesB_indSU3();

    // map to store all the parameters used in amplitudes
    map<string, vector<string>> channelParameters;

    // map with parameter values (filled from BAT each iteration)
    map<string, double> parameterValues;

    map<string, double> obs;

    double LogLikelihood(const vector<double> &parameters);

    void PrintHistogram();

    // channel lists
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

    // physical constants
    const double m_Bp = 5.27941;
    const double m_Bd = 5.27972;
    const double m_Bs = 5.36693;
    const double tau_Bp = 1.638e-12;
    const double tau_Bd = 1.517e-12;
    const double tau_Bs = 1.5195e-12;
    const double h_t = 6.582119569e-25;
    const double G_F = 1.1663788e-5;

    map<string, double> mesonMasses = {
        {"jpsi", 3.09690}, {"k0", 0.49761}, {"k0s", 0.49761}, {"k0l", 0.49761},
        {"k0b", 0.49761}, {"p0", 0.134976}, {"kp", 0.49367}, {"pp", 0.139570},
        {"om", 0.78266}, {"ph", 1.019461}, {"rh", 1.465},
        {"rho0", 0.77526}, {"rhop", 0.77511}, {"kst", 0.845}, {"kbst", 0.845},
        {"dp", 1.86966}, {"dm", 1.86966}, {"d0b", 1.86484}, {"d0", 1.86484},
        {"dsp", 1.96835}, {"dsm", 1.96835}, {"eta", 0.547862}, {"etap", 0.95778}};

    pair<string, pair<string, string>> parseChannel(const string &channel) const;
    double getBMesonLifetime(const string &bMeson) const;
    double getBMesonMass(const string &bMeson) const;

    void DefineParameters(const string &channel);
    map<string, double> DeclareParameters();

    // Direct-lookup getter (no delta-chain resolution)
    TComplex getPar(const string &name) const;

    void SetParameterValue(const string &paramName, double value);
    double getParameterValue(const string &paramName) const;

    map<string, dato> getMeas() const { return meas; }
    void setMeas(const map<string, dato> &newMeas) { meas = newMeas; }

    map<string, double> CalculatePolarizations(
        pair<TComplex, TComplex> &amplitude_0,
        pair<TComplex, TComplex> &amplitude_paral,
        pair<TComplex, TComplex> &amplitude_perp);

    void compute_decay_amplitudes(const string &channel);

    double CalculateBR(TComplex amplitude, TComplex amplitude_conj, const string &channel) const;
    double CalculateAcp(const TComplex &amplitude, const TComplex &conjugate_amplitude) const;
    double CalculateAlpha(const TComplex &amplitude, const TComplex &conjugate_amplitude, const string &channel);
    double CalculateC(const TComplex &amplitude, const TComplex &conjugate_amplitude, const string &channel);
    pair<double, double> CalculateS(const TComplex &amplitude, const TComplex &conjugate_amplitude, const string &channel);
    tuple<double, double, double> CalculatePhiAndLambda(const TComplex &amplitude, const TComplex &conjugate_amplitude, const string &channel);
    pair<vector<string>, string> extractChannelFromCorrKey(const string &corr_key);
    map<string, double> getPolarizationParams(const string &channel, const map<string, pair<TComplex, TComplex>> &amplitude_map);

    double Calculate_UncorrelatedObservables(map<string, pair<TComplex, TComplex>> &amplitude_map);
    double Calculate_CorrelatedObservables(map<string, pair<TComplex, TComplex>> &amplitude_map);

    // SU(3) and EWP penalty functions
    double CalculateSU3Penalty(double sigma) const;
    double CalculateEWPPenalty() const;

    void MCMCUserIterationInterface();
    void SaveHistograms(const string &filename);
    void PrintObservablePulls(const string &filename);

private:
    TComplex lam_bs_c, lam_bs_u, lam_bd_c, lam_bd_u;
    TComplex lamst_bs_c, lamst_bs_u, lamst_bd_c, lamst_bd_u;

    double ewp_limit = 0.0;
    double su3_sigma = 0.3;    // fixed value (used when su3_sigma_is_free = false)
    bool su3_sigma_is_free = false;

    // List of adjacent SU(3)-related amplitude base-name pairs
    // (each entry is a pair of base-names, getPar() is called on each)
    vector<pair<string, string>> su3Pairs;

    // List of EW-penguin parameter base names (for Gaussian prior)
    vector<string> ewpParamBaseNames;

    map<string, dato> meas;
    map<string, dato> newmeas;
    Histos histos;
    map<string, CorrelatedGaussianObservables> corrmeas;
    map<string, vector<string>> corrmeas_channels;
    unordered_set<string> referenceAmplitudes;
    map<string, pair<TComplex, TComplex>> amplitude_map;

    // Helper: add a SU(3) pair (base names, without _re/_im)
    // For vector channels (isVector=true) all three polarizations are registered
    void addSU3Pair(const string &param1, const string &param2, bool isVector = false)
    {
        if (!isVector) {
            su3Pairs.push_back({param1, param2});
        } else {
            for (const string &pol : {"_0", "_paral", "_perp"}) {
                su3Pairs.push_back({param1 + pol, param2 + pol});
            }
        }
    }

    // Helper: register an EW penguin base name (only once)
    void registerEWP(const string &baseName)
    {
        if (find(ewpParamBaseNames.begin(), ewpParamBaseNames.end(), baseName) == ewpParamBaseNames.end())
            ewpParamBaseNames.push_back(baseName);
    }

    string addPolarizationSuffix(string amplitude, string suffix) const
    {
        size_t pos_re = amplitude.find("_re");
        size_t pos_im = amplitude.find("_im");
        if (pos_re != string::npos)
            amplitude.insert(pos_re, suffix);
        else if (pos_im != string::npos)
            amplitude.insert(pos_im, suffix);
        else {
            cerr << "Target substring not found!" << endl;
            exit(1);
        }
        return amplitude;
    }

    int evaluatedevts;

    // Add a direct independent amplitude parameter (BAT variable)
    void addAmplitudeParameter(const string &amplitudeName, const double &lowerLimit, const double &upperLimit, bool isVectorChannel = false)
    {
        if (!isVectorChannel) {
            if (!referenceAmplitudes.count(amplitudeName)) {
                AddParameter(amplitudeName, lowerLimit, upperLimit);
                referenceAmplitudes.emplace(amplitudeName);
            }
        } else {
            vector<string> polarizations = {"_0", "_paral", "_perp"};
            for (const auto &pol : polarizations) {
                string fullName = addPolarizationSuffix(amplitudeName, pol);
                if (!referenceAmplitudes.count(fullName)) {
                    AddParameter(fullName, lowerLimit, upperLimit);
                    referenceAmplitudes.emplace(fullName);
                }
            }
        }
    }
};

#endif // GOLDENMODESB_INDSU3_H
