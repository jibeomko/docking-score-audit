# What a docking score can and cannot tell you

A computed binding score is asked to stand in for three different things: is this
a binder, how tightly does it bind, and will it work. This is a measurement of
how far it gets, on a compound set where all three answers are already known
from experiment.

> **Companion review** — Hwang Y, Kim J, **Ko JB**, Han Y, Kim Y, Lee KH, Jang S.
> *From affinity to kinetics: SPR as the analytical backbone of AI-driven drug discovery.*
> **TrAC Trends in Analytical Chemistry** 204:119089 (2026).
> [doi:10.1016/j.trac.2026.119089](https://doi.org/10.1016/j.trac.2026.119089) · CC BY-NC 4.0
>
> The review argues that affinity-centric pipelines discard the kinetic
> quantities that drive in vivo behaviour. This repository tests that claim
> quantitatively, on the review's own case study, with data the review did not
> contain.

## The short version

Ten adenosine A<sub>2A</sub> receptor agonists have published association and
dissociation rate constants, affinities and functional efficacies, all measured
in a single study. Rank them by each criterion a project team could actually
use, take the top three, and look at what you selected.

| selection criterion | compounds picked | mean efficacy |
|---|---|---|
| **residence time, from binding kinetics** | UK432097, CGS21680, LUF5549 | **102.0** |
| random pick | — | 77.6 |
| GNINA CNNaffinity | LUF5835, CGS21680, LUF5834 | 67.0 |
| Vina | CGS21680, LUF5834, LUF5835 | 67.0 |
| affinity, 1/K<sub>i</sub> | LUF5835, LUF5834, LUF5833 | 51.7 |

Selecting on the kinetic parameter returns the best possible set of three, 102.0
being the maximum achievable from this set. Selecting on measured affinity does
26 points worse than chance. Both docking scores do 11 points worse than chance.

The cause is visible in the raw table: the three compounds with the shortest
residence times score highest, and the compound with the longest residence time
and the highest efficacy scores near the bottom.

## The set

Ground truth is Guo et al., *Br J Pharmacol* 2012;166:1846-1859, the study the
companion review cites for the efficacy/residence-time relationship. Compound
identity was assigned by matching **both** rate constants against ChEMBL records
for ADORA2A. All ten matched k<sub>on</sub> and k<sub>off</sub> simultaneously,
a two-parameter agreement that fixes the mapping, so the structures docked here
are ChEMBL's rather than redrawn by hand.

Docking used GNINA v1.3.3 against three agonist-bound A<sub>2A</sub> structures
(3QAK, 2YDV, 2YDO), three seeds each, 90 prospective runs.

## Validation, before any score was interpreted

**The transcribed ground truth reproduces the published statistics.** This is the
check that the table was copied faithfully, run before the docking numbers were
looked at.

| relationship | reproduced | published |
|---|---|---|
| efficacy vs log residence time | r² = 0.90 | 0.90 |
| efficacy vs log K<sub>i</sub> | r² = 0.14 | 0.13 |
| cAMP efficacy vs log residence time | r² = 0.73 | 0.74 |
| cAMP efficacy vs log K<sub>i</sub> | r² = 0.11 | 0.10 |

**Redocking.** Each structure had its own co-crystal ligand docked back in. NECA
in 2YDV recovered to 0.23 Å and adenosine in 2YDO to 0.57 Å. UK-432097 in 3QAK
failed at 6.76 Å.

**Search failure or scoring failure?** The failed case was rerun at twice and
four times the sampling budget. The best recoverable pose stayed at exactly
6.48 Å at every budget while the top-scored pose got *worse*, and the score of
that wrong pose was the highest in the whole set. The crystal pose is not the
scoring optimum, so no additional search reaches it. Across the 90 prospective
runs the seed-to-seed score spread stayed under 0.35 log units for every
compound and under 0.10 for five of ten, so the sampling itself had converged.

## Result

| predictor | target | Spearman ρ | p | r² | 95 % CI on r² |
|---|---|---|---|---|---|
| CNNaffinity | pK<sub>i</sub> | +0.43 | 0.21 | 0.02 | 0.00–0.73 |
| CNNaffinity | log residence time | **−0.38** | 0.28 | 0.16 | 0.00–0.91 |
| CNNaffinity | efficacy | **−0.43** | 0.25 | 0.24 | 0.00–0.96 |
| Vina | pK<sub>i</sub> | −0.09 | 0.80 | 0.06 | 0.00–0.82 |
| Vina | log residence time | +0.19 | 0.60 | 0.42 | 0.01–0.88 |

Three things follow.

**It does not rank affinity.** Against experimental pK<sub>i</sub>, the quantity
the scoring function is trained to reproduce, r² is 0.02.

**It does not reach kinetics, and leans the wrong way.** Both kinetically
relevant correlations carry a negative sign. A higher computed affinity goes
with a *shorter* residence time and *lower* efficacy. For Vina the apparent
r² of 0.42 against residence time comes with ρ = +0.19 and p = 0.60, so it is
leverage from the extremes rather than a trend, and its sign is unhelpful too.

**The sample size is not the excuse.** Nothing reaches significance and every
interval spans most of the unit range, because n = 10. But the same ten
compounds support r² = 0.90 for efficacy against residence time at p < 0.0001.
The experiment resolves at this n. The score does not.

One structural caveat is worth stating plainly. The only receptor whose pocket
accommodates the 778 Da UK-432097 is 3QAK, the one that failed redocking, and
the two that passed cannot fit it without clashing. For the single most
interesting compound in the set there is no structure that is both validated and
usable.

Figures: `figures/a2a/FigA_experiment.pdf` (the published relationship),
`FigB_score.pdf` (the score against the same three axes),
`FigC_sampling.pdf` (sampling budget versus redocking RMSD).

## Why this sits next to an SPR review

The review's claim is that a single affinity number is the wrong target, because
k<sub>on</sub>, k<sub>off</sub> and residence time carry the decision-relevant
information. This audit reaches the same place from the computational side and
adds a sharper version of it. The score does not recover affinity, so the
affinity-centric pipeline is weaker than advertised even on its own terms. And
affinity was the wrong target anyway: on this set the kinetic parameter selects
the optimal compounds while affinity-based selection performs worse than chance.

A static score cannot express a rate constant. No rescoring of these pipelines
separates two compounds of equal computed affinity and seventy-fold different
residence time. Closing that gap takes a measurement, which is what the review
assigns to surface plasmon resonance.

## Layout

```
pipeline/a2a/   dataset, preparation, docking, analysis, figures
results/a2a/    derived tables and metrics only
figures/a2a/    publication-resolution PDF panels
```

Raw prediction output is not committed. The run is 102 docking jobs and takes
about 1.3 h of CPU wall time. Everything in `results/` is derived from it and is
sufficient to reproduce every number above.

## Reproducing

```bash
export A2A_WORK=~/a2a_work           # raw output goes here, outside the repo (default)
export GNINA=/path/to/gnina          # default: gnina on PATH
export OBABEL=/path/to/obabel        # default: obabel on PATH
export CUDNN_LIB=/path/to/cudnn/lib  # only if your gnina build needs it

python3 pipeline/a2a/dataset.py   # self-check of the transcribed ground truth
python3 pipeline/a2a/prep.py      # ChEMBL structures, receptors, 3D ligands
python3 pipeline/a2a/dock.py      # sampling sweep, redocking, 90 docking runs
python3 pipeline/a2a/analyze.py   # correlations and the selection simulation
python3 pipeline/a2a/figure.py
```

GNINA v1.3.3, RDKit, Open Babel, SciPy. Receptor structures from the PDB; ligand
structures and kinetic records from ChEMBL; kinetic and efficacy ground truth
from Guo et al. 2012.

## Scope and limits

Everything here is computed from existing experimental structures. No compound
was synthesised and no binding was measured. The set has n = 10, the size of the
published kinetic study, and the correlation estimates are correspondingly wide.
The conclusion drawn is a comparison against a relationship that *is* resolvable
at that n, not a claim that the score is significantly anticorrelated with
anything. Results are not evidence of activity for any compound and are not fit
for clinical or regulatory use.

## Licence

Code is MIT (`LICENSE`). The companion review is CC BY-NC 4.0 and is cited, not
redistributed. Ground-truth values are transcribed from the published literature
with the source named at the point of use.

## Author

Ji Beom Ko · co-author of the companion review.
