# What a docking score can and cannot tell you

A docking program takes a drug-like molecule and a protein, tries many ways of
fitting one into the other, and returns a single number. That number is meant to
say how well the molecule sticks. It is cheap to compute, which is why virtual
screening leans on it, and it ends up standing in for three different questions
at once: is this a binder, how tightly does it bind, and will it work.

Those are three questions, not one, and answering the first does not answer the
other two. This is a measurement of how far the score gets, on a compound set
where all three answers are already known from experiment. Terms that are
specific to the field are collected in [Terms](#terms).

> **Related review** — Hwang Y<sup>†</sup>, Kim J<sup>†</sup>, **Ko JB**<sup>†</sup>, Han Y, Kim Y, Lee KH<sup>\*</sup>, Jang S<sup>\*</sup>.
> *From affinity to kinetics: SPR as the analytical backbone of AI-driven drug discovery.*
> **TrAC Trends in Analytical Chemistry** 204:119089 (2026).
> [doi:10.1016/j.trac.2026.119089](https://doi.org/10.1016/j.trac.2026.119089) · CC BY-NC 4.0 ·
> <sup>†</sup> co-first authors, <sup>\*</sup> corresponding authors
>
> The review argues that affinity-centric pipelines discard the kinetic
> quantities that drive in vivo behaviour. This repository tests that claim
> quantitatively, on the A<sub>2A</sub> example the review uses (Section 2.2,
> Fig. 1B–D), with data the review did not contain.
>
> This analysis was done after the review was published. It is not part of the
> review and was not reviewed by its co-authors.

## The short version

Ten agonists of the adenosine A<sub>2A</sub> receptor, molecules that switch that
receptor on, were characterised in a single study. For each one it reported how
fast the compound binds, how fast it lets go, how tightly it holds, and how large
a response it produces. That last quantity, efficacy, is the one a drug is
ultimately judged on.

So rank the ten by each criterion a project team could actually use, take the top
three, and look at what you selected.

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

## Terms

| term | what it means |
|---|---|
| agonist | a molecule that switches a receptor on, rather than blocking it. |
| affinity, K<sub>i</sub> | how tightly a compound holds the receptor once bound. A lower K<sub>i</sub> is tighter. pK<sub>i</sub> is its negative logarithm, so a higher pK<sub>i</sub> is tighter. |
| k<sub>on</sub> | how fast the compound finds the receptor and binds. |
| k<sub>off</sub> | how fast it falls off again. |
| residence time | how long a single molecule stays bound, which is 1/k<sub>off</sub>. A compound can bind tightly and still let go quickly. |
| efficacy | how large a functional response the bound compound produces, as a percentage of a reference agonist. |
| docking | fitting a 3D ligand structure into a receptor pocket by computer and scoring the fit. The fitted arrangement is a pose. |
| GNINA, Vina | docking programs. The Vina score is a physics-style function; GNINA's CNNaffinity is a neural network trained on known structures and their measured affinities. |
| PDB, ChEMBL | public databases of experimental protein structures and of measured bioactivity. A code such as 3QAK names one structure. |
| redocking | taking a ligand out of the crystal structure it came from and docking it back in, as a positive control. |
| RMSD | how far the docked pose sits from the experimental one, in ångström. Under roughly 2 Å is normally called a success. |
| Spearman ρ, r², p | rank correlation from +1 to −1, where the sign says which way the trend runs; fraction of variance explained, from 0 to 1; and the probability of seeing a correlation this strong if there were none. |
| SPR | surface plasmon resonance, the experimental method that measures k<sub>on</sub> and k<sub>off</sub> directly. |

## The set

Ground truth, meaning the measured values the predictions get judged against, is
Guo et al., *Br J Pharmacol* 2012;166:1846-1859, the study the review cites for
the efficacy/residence-time relationship. Compound identity was assigned by
matching **both** rate constants against ChEMBL records for ADORA2A, the gene
that encodes the A<sub>2A</sub> receptor. All ten matched k<sub>on</sub> and
k<sub>off</sub> simultaneously, a two-parameter agreement that fixes the mapping,
so the structures docked here are ChEMBL's rather than redrawn by hand.

Docking used GNINA v1.3.3 against three agonist-bound A<sub>2A</sub> structures
(3QAK, 2YDV, 2YDO). Docking search is stochastic, so each structure was run from
three random seeds, giving 90 prospective runs.

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

**Redocking.** Each structure had the ligand it was crystallised with docked back
in, which is the check that the method can recover an answer already known. NECA
in 2YDV recovered to 0.23 Å and adenosine in 2YDO to 0.57 Å. UK-432097 in 3QAK
failed at 6.76 Å.

**Search failure or scoring failure?** A docking run can miss for two reasons: it
never samples the right pose, or it samples it and scores it badly. The failed
case was rerun at twice and four times the sampling budget. The best recoverable
pose stayed at exactly 6.48 Å at every budget while the top-scored pose got
*worse*, and the score of that wrong pose was the highest in the whole set. The
crystal pose is not the scoring optimum, so no additional search reaches it.
Across the 90 prospective runs the seed-to-seed score spread stayed under 0.35
log units for every compound and under 0.10 for five of ten, so the sampling
itself had converged.

## Result

| predictor | target | Spearman ρ | p | r² | 95 % CI on r² |
|---|---|---|---|---|---|
| CNNaffinity | pK<sub>i</sub> | +0.43 | 0.21 | 0.02 | 0.00–0.73 |
| CNNaffinity | log residence time | **−0.38** | 0.28 | 0.16 | 0.00–0.91 |
| CNNaffinity | efficacy | **−0.43** | 0.25 | 0.24 | 0.00–0.96 |
| Vina | pK<sub>i</sub> | −0.09 | 0.80 | 0.06 | 0.00–0.82 |
| Vina | log residence time | +0.19 | 0.60 | 0.42 | 0.01–0.88 |

Read the sign before the magnitude. A negative ρ means the score ranks the
compounds in the opposite order to the experiment, which is worse than ranking
them at random. Three things follow.

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

One structural caveat is worth stating plainly. UK-432097 is large for a
drug-like molecule at 778 Da, and the only receptor whose pocket accommodates it
is 3QAK, the one that failed redocking; the two that passed cannot fit it without
clashing. For the single most interesting compound in the set there is no
structure that is both validated and usable.

Figures: `figures/a2a/FigA_experiment.pdf` (the published relationship),
`FigB_score.pdf` (the score against the same three axes),
`FigC_sampling.pdf` (sampling budget versus redocking RMSD).

## Why this sits next to a surface plasmon resonance review

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
assigns to SPR.

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

Needs GNINA v1.3.3, RDKit, Open Babel and SciPy. The last three install from
conda-forge or pip; GNINA ships as a release binary or builds from source at
<https://github.com/gnina/gnina>.

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

Receptor structures come from the PDB; ligand structures and kinetic records from
ChEMBL; kinetic and efficacy ground truth from Guo et al. 2012.

## Scope and limits

Everything here is computed from existing experimental structures. No compound
was synthesised and no binding was measured. The set has n = 10, the size of the
published kinetic study, and the correlation estimates are correspondingly wide.
The conclusion drawn is a comparison against a relationship that *is* resolvable
at that n, not a claim that the score is significantly anticorrelated with
anything. Results are not evidence of activity for any compound and are not fit
for clinical or regulatory use.

## Licence

Code is MIT (`LICENSE`). The review is CC BY-NC 4.0 and is cited, not
redistributed. Ground-truth values are transcribed from the published literature
with the source named at the point of use.

## Author

Ji Beom Ko · co-author of the review.
