#include "LikelihoodWeights.h"
#include <TBranch.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <iomanip>
#include <sstream>

LikelihoodWeightOutput& LikelihoodWeightOutput::operator=(const LikelihoodWeightOutput& other)
{
    if (this != &other) {
        enabled = other.enabled; capturing = schemaChanged = false;
        entries.clear(); file = nullptr; tree = schema = nullptr;
        parameterNames.clear(); parameterValues.clear();
    }
    return *this;
}

std::string LikelihoodWeightOutput::SafeName(const std::string& raw)
{
    // Injective ASCII encoding. Escape literal "_x" too, so e.g. a-b and
    // a_x2d_b cannot silently collide after punctuation is encoded.
    const char* hex = "0123456789abcdef";
    std::string safe;
    for (size_t i = 0; i < raw.size(); ++i) {
        const unsigned char c = raw[i];
        const bool letter = (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
        const bool digit = c >= '0' && c <= '9';
        if (letter || (digit && i) || (c == '_' && (i+1 == raw.size() || raw[i+1] != 'x')))
            safe += c;
        else { safe += "_x"; safe += hex[c/16]; safe += hex[c%16]; safe += '_'; }
    }
    if (safe.empty()) throw std::invalid_argument("Empty likelihood-weight name");
    return safe;
}

void LikelihoodWeightOutput::Initialize(TFile* destination, const std::string& model,
                                        const std::vector<std::string>& parameters)
{
    if (!enabled) return;
    if (!destination || !destination->IsOpen() || !destination->IsWritable())
        throw std::runtime_error("Likelihood weights require writable MCMC ROOT output");
    const std::string name = model + "_Observables";
    if (destination == file && destination->Get(name.c_str()) == tree && tree) return;
    if (destination->Get(name.c_str()) || destination->Get((model+"_LikelihoodWeightInfo").c_str()))
        throw std::runtime_error("Likelihood-weight trees already exist in this output");
    file = destination;
    TDirectory::TContext context(file);
    tree = new TTree(name.c_str(), "Unnormalized log-likelihood terms at saved production states");
    tree->Branch("Chain", &chain, "chain/i");
    tree->Branch("Iteration", &iteration, "iteration/i");
    tree->Branch("Sweep", &sweep, "sweep/i");
    tree->Branch("Phase", &phase, "phase/I");
    tree->Branch("Replica", &replica, "replica/i");
    tree->Branch("LogLikelihood", &likelihood, "log_likelihood/D");
    tree->Branch("LogPrior", &prior, "log_prior/D");
    tree->Branch("LogProbability", &probability, "log_probability/D");
    tree->Branch("WeightSum", &weightSum, "weight_sum/D");
    parameterNames = parameters;
    parameterValues.assign(parameters.size(), 0.);
    for (unsigned p = 0; p < parameters.size(); ++p) {
        if (tree->GetBranch(parameters[p].c_str()))
            throw std::runtime_error("Sampling parameter collides with likelihood metadata: " + parameters[p]);
        tree->Branch(parameters[p].c_str(), &parameterValues[p], (parameters[p]+"/D").c_str());
    }
    tree->SetAutoSave(10000000); // Bound memory and make in-flight ROOT recovery useful.
    schema = new TTree((model+"_LikelihoodWeightInfo").c_str(), "Likelihood-weight branch definitions");
    schema->Branch("branch", &infoBranch); schema->Branch("name", &infoName);
    schema->Branch("kind", &infoKind); schema->Branch("group", &infoGroup);
    schema->Branch("definition", &infoDefinition);
    schema->Branch("additive", &infoAdditive, "additive/O");
    for (auto& item : entries) AddBranches(item.second);
}

void LikelihoodWeightOutput::BeginPoint()
{
    capturing = enabled;
    if (!capturing) return;
    for (auto& item : entries) { item.second.value = 0.; item.second.active = false; }
}

void LikelihoodWeightOutput::AddBranches(Entry& entry)
{
    if (!tree) return;
    if (tree->GetBranch(entry.branch.c_str()))
        throw std::runtime_error("Duplicate likelihood-weight branch: " + entry.branch);
    TDirectory::TContext context(file);
    auto value = tree->Branch(entry.branch.c_str(), &entry.value, (entry.branch+"/D").c_str());
    const std::string active = entry.branch + "_active";
    auto flag = tree->Branch(active.c_str(), &entry.active, (active+"/O").c_str());
    if (!value || !flag) throw std::runtime_error("Cannot create likelihood-weight branches");
    // A term can be conditionally absent (e.g. an undefined zero-amplitude ratio).
    // New branches are backfilled with inactive zeroes, never stale weights.
    const double saved = entry.value; const bool wasActive = entry.active;
    entry.value = 0.; entry.active = false;
    for (Long64_t i = 0; i < tree->GetEntries(); ++i)
        if (value->Fill() < 0 || flag->Fill() < 0) throw std::runtime_error("Cannot backfill likelihood-weight branch");
    entry.value = saved; entry.active = wasActive;
    infoBranch = entry.branch; infoName = entry.name; infoKind = entry.kind;
    infoGroup = entry.group; infoAdditive = entry.additive;
    infoDefinition = entry.additive ? "additive log-likelihood term" :
        "logL_block(r)-logL_block(r with this residual set to zero); not additive or marginal leave-one-out";
    if (schema->Fill() < 0) throw std::runtime_error("Cannot write likelihood-weight definitions");
    schemaChanged = true;
}

void LikelihoodWeightOutput::Record(const std::string& name, double value,
                                    const std::string& kind, const std::string& group)
{
    if (!capturing) return;
    if (!std::isfinite(value)) throw std::runtime_error("Nonfinite likelihood weight: " + name);
    const Key key(kind, group, name);
    auto found = entries.find(key);
    if (found == entries.end()) {
        Entry entry;
        entry.name = name; entry.kind = kind; entry.group = group;
        entry.additive = kind != "correlated_member";
        entry.branch = (entry.additive ? SafeName(name) : SafeName(group)+"__"+SafeName(name)) + "_weight";
        for (const auto& item : entries)
            if (item.second.branch == entry.branch)
                throw std::runtime_error("Ambiguous likelihood-weight identifiers: " + entry.branch);
        found = entries.emplace(key, entry).first;
        AddBranches(found->second);
    }
    found->second.value += value;
    found->second.active = true;
}

double LikelihoodWeightOutput::Gaussian(const std::string& name, double residual, double sigma)
{
    const double value = -0.5 * (residual * residual / (sigma * sigma));
    Record(name, value);
    return value;
}

void LikelihoodWeightOutput::Fill(unsigned c, unsigned i, unsigned s, unsigned r,
                                 double ll, double recomputed, double lp, double probabilityValue,
                                 const std::vector<double>& parameters)
{
    if (!enabled) return;
    if (!tree || capturing) throw std::logic_error("Likelihood output requires an initialized tree and a completed point capture");
    if (parameters.size() != parameterValues.size())
        throw std::runtime_error("Likelihood-output parameter dimensions differ from the MCMC state");
    // Never resize the buffers attached to ROOT branches.
    std::copy(parameters.begin(), parameters.end(), parameterValues.begin());
    weightSum = 0.;
    for (const auto& item : entries) if (item.second.additive) weightSum += item.second.value;
    const double tolerance = 1.e-10 * std::max(1., std::abs(ll));
    if (!std::isfinite(ll) || !std::isfinite(lp) || !std::isfinite(probabilityValue) ||
        !std::isfinite(weightSum) || !std::isfinite(recomputed) ||
        std::abs(weightSum-ll) > tolerance || std::abs(recomputed-ll) > tolerance)
    {
        std::ostringstream error;
        error << std::setprecision(17) << "Saved likelihood contributions do not match the current MCMC state"
              << " (chain=" << c << ", iteration=" << i << ", sweep=" << s
              << "): sum=" << weightSum << ", recomputed=" << recomputed << ", cached=" << ll
              << ", prior=" << lp << ", probability=" << probabilityValue;
        throw std::runtime_error(error.str());
    }
    chain = c; iteration = i; sweep = s; replica = r;
    likelihood = ll; prior = lp; probability = probabilityValue;
    if (schemaChanged) {
        if (schema->AutoSave("SaveSelf") < 0) throw std::runtime_error("Cannot save likelihood-weight definitions");
        schemaChanged = false;
    }
    if (tree->Fill() < 0) throw std::runtime_error("Cannot write likelihood weights");
}

void GoldenModesLikelihoodModel::SetWriteMCMCWeights(bool flag)
{
    if (GetPhase() != kUnsetPhase || !fMCMCStates.empty() || fMCMCOutputFile)
        throw std::logic_error("Configure likelihood-weight output before initializing MCMC output");
    likelihoodWeights.Enable(flag);
}

void GoldenModesLikelihoodModel::InitializeMarkovChainTree(bool replacetree, bool replacefile)
{
    BCModel::InitializeMarkovChainTree(replacetree, replacefile);
    std::vector<std::string> parameters;
    if (likelihoodWeights.Enabled())
        for (unsigned p = 0; p < GetNParameters(); ++p) parameters.push_back(GetParameter(p).GetSafeName());
    likelihoodWeights.Initialize(fMCMCOutputFile, GetSafeName(), parameters);
}

double GoldenModesLikelihoodModel::EvaluateOutputLikelihood(const std::vector<double>& parameters)
{
    LikelihoodWeightOutput::Capture capture(likelihoodWeights);
    return LogLikelihood(parameters);
}

double GoldenModesLikelihoodModel::CorrelatedWeight(const std::string& group,
    const CorrelatedGaussianObservables& measurement, const TVectorD& prediction,
    const std::vector<std::string>& names)
{
    if (!likelihoodWeights.Capturing()) return measurement.logweight(prediction);
    if (names.size() != measurement.getNObs())
        throw std::runtime_error("Correlated likelihood-weight names have the wrong size: " + group);
    const auto weights = measurement.computeLeaveOneOutWeights(prediction);
    likelihoodWeights.Record(group, weights[0], "correlated_block");
    for (unsigned i = 0; i < names.size(); ++i)
        likelihoodWeights.Record(names[i], weights[i+1], "correlated_member", group);
    return weights[0];
}

void GoldenModesLikelihoodModel::WriteLikelihoodWeights(unsigned chain, double recomputed)
{
    if (!likelihoodWeights.Enabled() || fMCMCPhase != kMainRun || !fMCMCFlagWriteChainToFile) return;
    const auto& state = fMCMCStates.at(chain);
    const auto& ids = GetReplicaIDs();
    likelihoodWeights.Fill(chain, state.iteration, fMCMCCurrentIteration,
        ids.empty() ? chain : ids.at(chain), state.log_likelihood, recomputed, state.log_prior, state.log_probability, state.parameters);
}
