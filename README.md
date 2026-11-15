# EvolveX for NOMINEE

<!-- Publication note: replace this comment with
> **Used in:** [Paper title](DOI or publication URL) — authors, year.
when the citation is available. -->

**Structure-guided antibody and nanobody sequence optimisation with explicit control of stability and mutation burden.**

This repository contains the modified [EvolveX](https://github.com/SwitchLab-VIB/EvolveX) implementation used in the optimisation stage of [NOMINEE](https://github.com/CubeJerry/dev). Starting from antibody–antigen complex structures and a defined set of mutable positions, EvolveX explores sequence variants using FoldX energy estimates.

This fork extends the original binding-centred search with a combined objective: favour predicted binding improvements, reward improved antibody stability, and discourage unnecessary divergence from the original sequence. Structural checks, mutation limits, and composition-dependent proposal weights provide additional constraints during the search.

[Overview](#overview) · [Scoring and acceptance](#scoring-and-acceptance) · [Design constraints](#design-constraints) · [Outputs and interpretation](#outputs-and-interpretation) · [Code and integration](#code-and-integration) · [Attribution](#attribution)

## Overview

EvolveX combines local Monte Carlo mutation steps with recombination between models derived from the same starting complex. FoldX evaluates proposed structures, and the search decides whether to retain each proposal or keep the current model.

The changes described here apply to the **genetic algorithm / Monte Carlo search** selected by `search_algorithm: GA`. The separate systematic search enumerates mutation combinations and does not use this acceptance rule.

| Component | Original EvolveX | This fork |
|---|---|---|
| Acceptance objective | Binding-energy changes | Combined binding, antibody stability, and mutation-count changes |
| Structural checks | Hard stability and intraclash filters | Revised stepwise filters and a stability limit relative to the original antibody |
| Sequence divergence | No explicit mutation-burden term in acceptance | A mutation-count penalty, a soft cap, and a hard cap |
| Mutation proposals | Fixed amino-acid proposal weights | Weights adjusted according to the current sequence composition |
| Search records | Structural and energetic properties of retained models | Additional proposed and retained objective changes, plus mutation burden |

## Scoring and acceptance

### Reference states and sign convention

Two reference states have different roles:

- **Current model:** the model being modified at a particular search step. Energy changes compare a proposal with its FoldX-generated wild-type counterpart for that step.
- **Original sequence:** the antibody sequence in the input complex, before search initialisation. Mutation burden and cumulative stability drift are measured against this reference.

For an energy-like score $G$, a proposed change is written as

```math
\Delta G = G_{\mathrm{proposed}} - G_{\mathrm{current}}.
```

Lower values are better: a negative binding change indicates improved predicted interaction energy, and a negative antibody stability change indicates improved predicted stability.

### Combined objective

The fork uses the following change in score for each proposal:

```math
\Delta E_{\mathrm{rewarded}} =
\Delta G_{\mathrm{bind,obj}}
+
0.25\,\Delta G_{\mathrm{stab}}
+
0.5\,\Delta n_{\mathrm{mut}}.
```

| Term | Meaning |
|---|---|
| $\Delta G_{\mathrm{bind,obj}}$ | Change in the binding objective, defined below |
| $\Delta G_{\mathrm{stab}}$ | Change in the FoldX stability energy of the antibody chain group |
| $\Delta n_{\mathrm{mut}}$ | Change in the number of mutable positions that differ from the original sequence |

The coefficients are fixed in the current implementation. The mutation-count coefficient converts a count into a contribution on the search-score scale; the combined objective is a heuristic score, not a measured thermodynamic free energy.

The mutation term uses a **count**, rather than an unscaled fraction:

```math
\Delta n_{\mathrm{mut}} =
N_{\mathrm{mutable}}
\left(f_{\mathrm{proposed}}-f_{\mathrm{current}}\right),
```

where $f$ is the fraction of mutable positions differing from the original sequence. Introducing one additional difference contributes $+0.5$ to the objective; restoring one original residue contributes $-0.5$. Replacing one non-original residue with another leaves the count unchanged.

This lets stability improvements or a reduction in mutation burden offset a modest worsening in predicted binding energy.

### Binding objective

With the standard FoldX calculation, the binding objective is simply the binding-energy change. When the optional water-aware calculation is enabled, the two changes are averaged:

```math
\Delta G_{\mathrm{bind,obj}} =
\begin{cases}
\Delta G_{\mathrm{bind}}, & \text{standard calculation}, \\[4pt]
\dfrac{\Delta G_{\mathrm{bind}}+\Delta G_{\mathrm{bind,water}}}{2},
& \text{water-aware calculation enabled}.
\end{cases}
```

### Metropolis acceptance

After the hard checks described below, the combined score enters a Metropolis acceptance step:

```math
P_{\mathrm{Metropolis}} =
\begin{cases}
1, & \Delta E_{\mathrm{rewarded}}\leq 0, \\[4pt]
\exp\!\left(-\dfrac{\Delta E_{\mathrm{rewarded}}}{0.5919}\right),
& \Delta E_{\mathrm{rewarded}}>0.
\end{cases}
```

Favourable moves pass this step, while unfavourable moves can still be retained probabilistically. The scale factor $0.5919$ is inherited from the original implementation.

**Passing the Metropolis step does not guarantee final acceptance:** the soft mutation cap can still reject a move that increases divergence.

## Design constraints

### Structural checks

The GA search applies structural checks from its first search iteration. A proposal is rejected if it:

- worsens antibody stability by more than $3.0$ in a single step;
- increases the antibody intraclash score by more than $3.0$ in a single step; or
- exceeds the following antibody stability limit relative to the original input:

```math
G_{\mathrm{stab,proposed}}
>
\max\!\left(
1.3\,G_{\mathrm{stab,original}},
G_{\mathrm{stab,original}}+10.0
\right).
```

The last check limits accumulated stability drift across many individually permissible steps. It uses the more permissive of the two bounds shown above. Intraclash is a FoldX score and should be distinguished from the energy terms.

### Mutation burden

Mutation burden is defined over the **mutable positions tracked by the search**, not over the full antibody sequence:

```math
f =
\frac{
\text{mutable positions differing from the original sequence}
}{
N_{\mathrm{mutable}}
}.
```

| Constraint | Behaviour |
|---|---|
| Soft cap: 25% | Above this fraction, an otherwise accepted move that further increases divergence faces an additional probabilistic rejection. |
| Hard cap: 50% | Proposals above this fraction are rejected. |

For divergence-increasing moves between the soft and hard caps, the probability of surviving the additional soft-cap check decreases linearly:

```math
P_{\mathrm{soft}} =
1-\frac{f_{\mathrm{proposed}}-0.25}{0.50-0.25}.
```

This additional check applies only when $f_{\mathrm{proposed}}>f_{\mathrm{current}}$ and $f_{\mathrm{proposed}}>0.25$. Moves that preserve or reduce divergence are exempt from the soft-cap check, but remain subject to the hard checks and Metropolis decision. An increase to exactly 50% has zero soft-cap survival probability.

### Composition-dependent proposals

The fork also adjusts how mutations are **proposed**, before they are evaluated. It monitors the combined fraction of **F, H, I, R, W, and Y** among the mutable positions and scales their baseline proposal weights as follows:

| Current fraction of these residues | Weight multiplier |
|---|---|
| Up to 15% | 1.0 |
| Above 15%, up to 30% | 0.7 |
| Above 30%, up to 40% | 0.5 |
| Above 40% | 0.3 |

The removed weight is redistributed proportionally among the other amino acids. Proposals still respect the allowed amino acids at each position.

This is a composition-control heuristic: it reduces further enrichment of the selected residue group without banning those residues. It is not an explicit prediction of aggregation, solubility, or experimental developability. Together with the state-dependent proposals and additional rejection rules, the search is an optimisation heuristic rather than an exact equilibrium sampler.

## Outputs and interpretation

For the GA search, the main outputs are:

| Output | Contents |
|---|---|
| `generated_models_info.csv` | Search-step records, retained model properties, mutation acceptance, and objective diagnostics |
| `model_PDB_files/` | Compressed PDB snapshots of the retained models |
| `EvolveX/` | Per-backbone working files and mutation-screening results |

The CSV distinguishes **proposed** score changes from **retained** changes. Columns prefixed with `proposed_` record the evaluated proposal even when it is rejected; the corresponding retained objective-delta columns are zero for rejected proposals. Mutation fraction and structural properties describe the retained state.

These scores support candidate generation and comparison. They are computational estimates: improved objective values alone do not establish improved experimental affinity, stability, or developability. Downstream assessment and experimental results should be reported separately.

## Code and integration

The package requires **Python 3.9 or later** and an external FoldX installation. Python dependencies and the `evolvex` command-line entry point are declared in [setup.py](setup.py).

For the integrated workflow, see [NOMINEE](https://github.com/CubeJerry/dev). For the original EvolveX usage documentation and examples, see the [upstream repository](https://github.com/SwitchLab-VIB/EvolveX). The upstream example configuration and example dataset files are not bundled in this fork.

## Attribution

EvolveX was developed by **Gabriel Cia and Rob Van Der Kant at SwitchLab**. This repository is a modified fork of [SwitchLab-VIB/EvolveX](https://github.com/SwitchLab-VIB/EvolveX); the upstream project provides the underlying FoldX-based design framework.

When reporting work that uses this implementation, distinguish the original EvolveX method from the modifications described here and identify the exact fork revision used. The associated NOMINEE publication citation can be added at the top of this README when available.

Code is distributed under the [Apache License 2.0](LICENSE). FoldX is obtained separately under its own licensing terms.
