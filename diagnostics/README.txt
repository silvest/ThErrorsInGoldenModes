GoldenModes-specific diagnostic definitions
==========================================

These configurations live with the model whose parameter names and
prior-region cuts they describe. Generic analysis code and statistical
regression tests remain in the bat_MPI repository.

Files:
  goldenmodes_partitions.json
    Separate binary views of the 14 parameter cuts used by the old region
    flags. Missing channel-specific parameters are optional and skipped.
    Probabilities from different partitions overlap and must not be added.

  goldenmodes_joint_regions.example.json
    Example disjoint BJPSIP bins combining the signs of two G2t amplitudes
    and E2t magnitude. Mixed sign configurations are reported as other.
    Adapt the bins to the physical modes of interest; these example bins
    do not automatically identify all posterior modes.

From the ThErrorsInGoldenModes repository, assuming bat_MPI is its sibling:

  python ../bat_MPI/diagnostics/analyze_mcmc.py results/run_cold_chains.root \
    --regions diagnostics/goldenmodes_partitions.json \
    --output results/run_diagnostics

For joint bins, replace --regions with:
  --regions diagnostics/goldenmodes_joint_regions.example.json

Install Python dependencies from ../bat_MPI/diagnostics/requirements.txt.
See ../PT_DIAGNOSTICS.txt (in this repository's root) for complete usage,
output descriptions and scientific interpretation.

Runtime region labels are defined in GoldenModesPTRegions.h in this directory.
Tests verify that those cuts match goldenmodes_partitions.json. See the parent
PT_DIAGNOSTICS.txt for streaming transport analysis and per-temperature output.
