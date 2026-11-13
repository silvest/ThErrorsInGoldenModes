#ifndef GOLDENMODES_LIKELIHOOD_WEIGHTS_H
#define GOLDENMODES_LIKELIHOOD_WEIGHTS_H

#include <BAT/BCModel.h>
#include "CorrelatedGaussianObservables.h"
#include <map>
#include <string>
#include <tuple>
#include <vector>

class TFile;
class TTree;

/// Opt-in output for saved current states. ROOT owns the trees; buffers stay here.
class LikelihoodWeightOutput {
public:
    struct Entry {
        std::string name, group, kind, branch;
        double value = 0.;
        bool active = false, additive = true;
    };
    using Key = std::tuple<std::string, std::string, std::string>;
    LikelihoodWeightOutput() = default;
    // Copy configuration, never ROOT objects or branch addresses.
    LikelihoodWeightOutput(const LikelihoodWeightOutput& other) : enabled(other.enabled) {}
    LikelihoodWeightOutput& operator=(const LikelihoodWeightOutput& other);
    void Enable(bool flag) { enabled = flag; }
    bool Enabled() const { return enabled; }
    bool Capturing() const { return capturing; }
    void Initialize(TFile* file, const std::string& model, const std::vector<std::string>& parameters = {});
    void BeginPoint();
    void EndPoint() { capturing = false; }
    void Record(const std::string& name, double value, const std::string& kind = "uncorrelated",
                const std::string& group = "");
    double Gaussian(const std::string& name, double residual, double sigma);
    void Fill(unsigned chain, unsigned iteration, unsigned sweep, unsigned replica,
              double likelihood, double recomputed, double prior, double probability,
              const std::vector<double>& parameters = {});
    const std::map<Key, Entry>& Entries() const { return entries; }
    static std::string SafeName(const std::string& raw);
    class Capture {
    public:
        explicit Capture(LikelihoodWeightOutput& output) : output(output) { output.BeginPoint(); }
        ~Capture() { output.EndPoint(); }
        Capture(const Capture&) = delete;
        Capture& operator=(const Capture&) = delete;
    private:
        LikelihoodWeightOutput& output;
    };
private:
    bool enabled = false, capturing = false, schemaChanged = false;
    std::map<Key, Entry> entries;
    std::vector<std::string> parameterNames;
    std::vector<double> parameterValues;
    TFile* file = nullptr;
    TTree *tree = nullptr, *schema = nullptr;
    unsigned chain = 0, iteration = 0, sweep = 0, replica = 0;
    int phase = 1;
    double likelihood = 0., prior = 0., probability = 0., weightSum = 0.;
    std::string infoBranch, infoName, infoKind, infoGroup, infoDefinition;
    bool infoAdditive = true;
    void AddBranches(Entry& entry);
};

/// Model-local facility; BAT's sampling and MPI worker protocol remain unchanged.
class GoldenModesLikelihoodModel : public BCModel {
public:
    void SetWriteMCMCWeights(bool flag = true);
    bool GetWriteMCMCWeights() const { return likelihoodWeights.Enabled(); }
    void InitializeMarkovChainTree(bool replacetree = false, bool replacefile = false) override;
protected:
    LikelihoodWeightOutput likelihoodWeights;
    double EvaluateOutputLikelihood(const std::vector<double>& parameters);
    double CorrelatedWeight(const std::string& group, const CorrelatedGaussianObservables& measurement,
                            const TVectorD& prediction, const std::vector<std::string>& names);
    void WriteLikelihoodWeights(unsigned chain, double recomputed);
};
#endif
