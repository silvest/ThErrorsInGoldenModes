#ifndef GOLDENMODES_PT_REGIONS_H
#define GOLDENMODES_PT_REGIONS_H
#include <BAT/BCEngineMCMC.h>
#include <BAT/BCParameter.h>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

// Model-owned binary cuts. The regression test compares these definitions with
// goldenmodes_partitions.json, used by the offline cold-chain analyzer.
namespace GoldenModesPT {
struct Cut { const char* parameter; char transform; double threshold; const char* yes; const char* no; };
inline const std::vector<Cut>& cuts() {
    static const std::vector<Cut> result = {
    {"G2t_scd_BJPSIP_arg", 'c', 0.0, "positive", "negative"},
    {"G2t_dcd_BJPSIP_arg", 'c', 0.0, "positive", "negative"},
    {"E2t_ccdd_BJPSIP_abs", 'i', 0.75, "small", "large"},
    {"EA1_sdcd_BPJPSI_arg", 'c', 0.0, "positive", "negative"},
    {"EA1_ddcd_BPJPSI_arg", 'c', 0.0, "positive", "negative"},
    {"EA1_sdcd_BPJPSI_abs", 'i', 1.5, "small", "large"},
    {"EA1_ddcd_BPJPSI_abs", 'i', 1.5, "small", "large"},
    {"EA2t_ccds_BJPSIP_arg", 's', 0.0, "positive", "negative"},
    {"E2t_ccds_BJPSIV_paral_arg", 'c', 0.0, "positive", "negative"},
    {"G2t_dcs_BJPSIV_0_arg", 'c', 0.0, "positive", "negative"},
    {"G2t_dcs_BJPSIV_paral_arg", 'c', 0.0, "positive", "negative"},
    {"G2t_dcs_BJPSIV_perp_arg", 'c', 0.0, "positive", "negative"},
    {"G1t_scd_BDDb_arg", 'c', 0.0, "positive", "negative"},
    {"A1_scdd_BDDb_arg", 'c', 0.0, "positive", "negative"},
    };
    return result;
}
inline std::vector<BCEngineMCMC::PTDiagnosticPartition> definitions(
        const BCEngineMCMC& model, std::vector<std::pair<unsigned, unsigned>>& indices)
{
    indices.clear();
    std::vector<BCEngineMCMC::PTDiagnosticPartition> result;
    for (unsigned c = 0; c < cuts().size(); ++c)
        for (unsigned i = 0; i < model.GetNParameters(); ++i)
            if (model.GetParameter(i).GetName() == cuts()[c].parameter) {
                indices.push_back({i, c});
                result.push_back({cuts()[c].parameter, {cuts()[c].yes, cuts()[c].no}});
                break;
            }
    return result;
}
inline void classify(const std::vector<std::pair<unsigned, unsigned>>& indices,
                     const std::vector<double>& parameters, std::vector<unsigned>& labels)
{
    labels.resize(indices.size());
    for (unsigned p = 0; p < indices.size(); ++p) {
        const Cut& cut = cuts()[indices[p].second];
        const double raw = parameters.at(indices[p].first);
        const double value = cut.transform == 'c' ? std::cos(raw) : cut.transform == 's' ? std::sin(raw) : raw;
        if (!std::isfinite(value)) throw std::runtime_error("Nonfinite GoldenModes PT region parameter");
        const bool yes = cut.transform == 'i' ? value < cut.threshold : value >= cut.threshold;
        labels[p] = yes ? 0 : 1;
    }
}
}
#endif
