# Per-observable likelihood weights

Add `--write-mcmc-weights` to a normal `mainGM` command. The option works with
both `goldenmodesB_indSU3` (`--indSU3`) and `goldenmodesB`, with ordinary MPI
Metropolis and with PT. It defaults to off. For example, add it to the existing
physics, warmup and temperature options:

```sh
mpiexec -n 3 ./mainGM --indSU3 --flagBJPSIP --minimal --Univariate \
  --pt-full-priors --pt-temperatures 1,2,4 --nChains 2 \
  --write-mcmc-weights --outfile fit_
```

PT output goes into `<outfile>cold_chains.root`; ordinary runs with this flag
also save parameter chains, in `<outfile>chains.root`. Rebuild `mainGM` after
updating the sources. No BAT or HEPfit change is needed for this feature.

## ROOT format and alignment

In the same file as the existing `<model>_mcmc` parameter tree, two trees are
added (the default model name is `model`):

- `<model>_Observables`: one row for every saved production state and cold
  chain, in exactly the same order as `<model>_mcmc`. The existing lag is applied.
  All ordinary chains are at the target temperature; PT saves the cold slots.
- `<model>_LikelihoodWeightInfo`: one definition per weight branch, with
  `branch`, `name`, `kind`, `group`, `definition`, and `additive` fields.

The weight tree includes `Chain`, `Iteration`, `Phase`, `Replica`, `Sweep`,
`LogLikelihood`, `LogPrior`, `LogProbability` and `WeightSum`. The first four
identifiers match the parameter tree row for row. `Sweep` is the production
sweep number; `Iteration` follows BAT's chain-state counter. `Replica` tracks
configurations through PT swaps, while `Chain` identifies the cold slot.
It also includes every sampling parameter, including fixed parameters, under
the same branch names as `<model>_mcmc`. This allows scripts that read the
parameter and weights from a single tree to work directly.

Weights are evaluated at the saved **current state**, after any swaps, including
repeated states from rejected proposals. Prerun, fixed settling and hot slots
are excluded. GoldenModes already reevaluates each saved cold state to fill
its histograms; the capture reuses that call. It does not add likelihood
evaluations to the sampler or change the MPI worker protocol, proposals or RNGs.
ROOT streams the rows rather than retaining a full chain in memory.

## Definitions

All values are **log weights**, with the same omitted normalization constants
as the model's likelihood. They are not exponentiated probabilities or importance
weights. Positive correlated-member contrasts are possible.

| `kind` | Branch convention | Meaning | Additive |
| --- | --- | --- | --- |
| `uncorrelated` | `<measurement-key>_weight` | Actual individual Gaussian log-likelihood term, including any phase wrapping used by the model | yes |
| `correlated_block` | `<dataset-key>_weight` | Full multivariate Gaussian log weight for that experimental dataset | yes |
| `correlated_member` | `<dataset-key>__<observable-name>_weight` | HEPfit-style residual-removal contrast within that dataset | no |
| `constraint` | `SU3Breaking_weight` | Full SU(3)-breaking penalty included in `goldenmodesB_indSU3::LogLikelihood` | yes |

Group-prefixed member names distinguish experiments that measure the same
observable. The schema retains the original measurement/dataset names.
Branch names preserve ordinary identifiers; punctuation, non-ASCII bytes,
initial digits and literal `_x` sequences are encoded as `_xHH_` as needed to
avoid ROOT identifier collisions. Use the schema's `branch` field rather than
reconstructing names. Each weight has a matching `<branch>_active` boolean.

The exporter records the actual terms added by the selected model and flags.
If a conditional term is not used at a state (e.g. a ratio skipped at zero
denominator), its weight is zero and its active flag is false. A term first
encountered later is backfilled with inactive zeros in earlier rows. Terms
never used by the likelihood at saved states have no branch. If a measurement
key is added more than once by the model, its recorded weight is their sum.

The additive terms, including the SU(3) constraint, sum to the saved
`LogLikelihood`. This is checked at every output row, together with agreement
between the cached state likelihood and its reevaluation. An inconsistent
capture stops the run instead of silently writing mismatched weights.

BAT priors, including CKM, mixing-angle and EWP Gaussian priors, are reflected
in `LogPrior`; they are not experimental likelihood terms. The existing EWP
prior is not added a second time as a likelihood penalty. `SU3Breaking_weight`
is explicitly labelled as a theoretical constraint so experimental information
can be studied separately. `LogProbability = LogLikelihood + LogPrior` follows
the untempered BAT convention, also for PT output.

## Correlated-member convention

This follows HEPfit's `computeLeaveOneOutWeights` implementation. Let
`r = prediction - data`, and let `P` be the inverse covariance matrix:

```text
block       = -1/2 r^T P r
member_i    = block(r) - block(r with r_i = 0)
            = -1/2 P_ii r_i^2 - r_i sum_{j != i} P_ij r_j
```

The other residuals and the full precision matrix are held fixed. This is
**not** the marginal likelihood obtained by deleting a measurement and
re-inverting its reduced covariance matrix. Member contrasts generally do not
sum to the block weight, because the cross terms occur in multiple contrasts.
For additive information accounting, use the full correlated block once.
Member contrasts support the same convention as the existing HEPfit output;
do not add them to the block or to `WeightSum`.

## Reading the output

For example, with Python and `uproot`:

```python
import numpy as np
import uproot

with uproot.open("fit_cold_chains.root") as root:
    parameters = root["model_mcmc"].arrays(library="np")
    weights = root["model_Observables"].arrays(library="np")
    schema = root["model_LikelihoodWeightInfo"].arrays(library="np")

for key in ("Chain", "Iteration", "Phase", "Replica"):
    np.testing.assert_array_equal(parameters[key], weights[key])
additive = [name for name, use in zip(schema["branch"], schema["additive"]) if use]
np.testing.assert_allclose(sum(weights[name] for name in additive),
                           weights["LogLikelihood"], rtol=1e-10, atol=1e-10)
experimental = [name for name, kind in zip(schema["branch"], schema["kind"])
                if kind in ("uncorrelated", "correlated_block")]
```

These arrays can feed a local-curvature or regression-based information proxy.
Log-weight values alone are not an exact Fisher matrix: the interpretation also
depends on how derivatives/curvature are estimated, parameter correlations,
priors and multimodality. This change exports the inputs; it does not impose
a particular Fisher-information estimator.

The existing HEPfit script can be run unchanged in its individual-observable
mode (note the capital `F` in the mounted directory):

```sh
python /Users/silvest/farm/SMEFT_Fits/psis_loo_weights.py \
  fit_cold_chains.root G2t_scd_BJPSIP_arg --top 20
```

Its FI proxy uses the stored per-state values directly. Its IS-LOO interpretation
applies to an additive individual likelihood or a full correlated block. For a
correlated-member branch it instead removes the residual-removal contrast, not
the measurement's marginal likelihood. The script does not currently consult
the `additive` or `kind` metadata: interpret those member rows as contrasts.
`SU3Breaking` is constraint sensitivity, not an experimental observable.
The script's source-location/category mode expects HEPfit `.conf` files;
GoldenModes defines its measurements in C++, so those source categories are
not automatically available. No original HEPfit script is modified here.

## Verification

With `GOLDENMODES_BUILD_TESTS=ON`, `likelihood_weights_ordinary` and
`likelihood_weights_pt` compare on/off trajectories using identical explicit
initial states and seeds. They test non-flat priors, analytic correlated blocks
and member contrasts, rejected states, replica exchanges, lag, discarded warmup,
metadata alignment, dynamically absent terms and rejection of inconsistent sums.
The number of likelihood evaluations is unchanged by enabling export.
`likelihood_weight_models` checks both physics models and all three SU(3)
penalty modes. `likelihood_weights_cli` checks real ordinary and PT driver runs.

```sh
ctest --test-dir build-cmake -R '^likelihood_weight' --output-on-failure
```

## Experimental-average initialization fix

`PDGAverage` now defaults to arithmetic averaging. Its angular flag was previously
uninitialized, so MPI ranks could construct different experimental averages.
Explicit circular averaging remains available. Rebuild for new runs; existing
chains retain their original likelihood evaluations. The `pdg_average` regression
test covers this initialization and explicit circular averaging.
