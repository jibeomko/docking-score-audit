#!/usr/bin/env python
"""A2A adenosine receptor agonist kinetic reference set.

Ground truth: Guo D, Mulder-Krieger T, IJzerman AP, Heitman LH.
"Functional efficacy of adenosine A2A receptor agonists positively correlated
to receptor residence time." Br J Pharmacol 2012;166:1846-1859. Table 3 for
binding kinetics and affinity, Results text for functional efficacy.

Compound identity was assigned by matching BOTH rate constants against ChEMBL
target CHEMBL251 (ADORA2A) kinetic records. All ten compounds matched kon and
koff simultaneously, so the ChEMBL ID -> compound name mapping is exact.
SMILES are taken from ChEMBL, not drawn by hand.

Units: kon M-1 min-1, koff min-1, RT min, KD and Ki nM, efficacy % of CGS21680.
"""

# name: (chembl_id, kon, koff, RT, kinetic_KD, Ki, eff_impedance, eff_cAMP)
COMPOUNDS = {
    "CGS21680": ("CHEMBL331372",  5.0e4, 0.020,  53.0, 380.0, 376.0, 100.0, 100.0),
    "NECA":     ("CHEMBL464859",  5.0e5, 0.030,  35.0,  58.0,  64.0,  None,  None),
    "UK432097": ("CHEMBL1096896", 5.0e5, 0.004, 250.0,   8.0,  22.0, 114.0, 115.0),
    "LUF5448":  ("CHEMBL3932873", 2.8e5, 0.060,  16.0, 225.0, 219.0,  83.0,  84.0),
    "LUF5549":  ("CHEMBL3981215", 2.4e6, 0.040,  24.0,  17.0,  24.0,  92.0,  71.0),
    "LUF5550":  ("CHEMBL3978846", 8.0e5, 0.090,  12.0, 110.0, 126.0,  63.0,  39.0),
    "LUF5631":  ("CHEMBL3983152", 8.0e5, 0.050,  21.0,  60.0,  44.0,  91.0,  67.0),
    "LUF5833":  ("CHEMBL124345",  8.5e6, 0.160,   6.3,  19.0,  17.0,  54.0,  38.0),
    "LUF5834":  ("CHEMBL122622",  1.1e7, 0.230,   4.2,  21.0,  16.0,  54.0,  50.0),
    "LUF5835":  ("CHEMBL122806",  1.6e7, 0.290,   3.4,  18.0,  15.0,  47.0,  58.0),
}

FIELDS = ("chembl_id", "kon", "koff", "RT", "kinetic_KD", "Ki",
          "eff_impedance", "eff_cAMP")

# Agonist-bound A2A structures used as docking receptors. All three carry a
# co-crystallised agonist, so each doubles as a redocking validation case.
RECEPTORS = {
    "3QAK": {"ligand": "UKA", "resolution": 2.71, "native": "UK432097",
             "note": "UK-432097 bound, agonist conformation"},
    "2YDV": {"ligand": "NEC", "resolution": 2.60, "native": "NECA",
             "note": "NECA bound, thermostabilised"},
    "2YDO": {"ligand": "ADN", "resolution": 3.00, "native": None,
             "note": "adenosine bound; adenosine is not in the kinetic set"},
}

# Correlations reported by Guo et al. 2012, reproduced here as the benchmark
# the docking scores are measured against.
PUBLISHED = {
    "eff_impedance_vs_logRT": {"r2": 0.90, "p": "<0.0001"},
    "eff_impedance_vs_logKi": {"r2": 0.13, "p": "0.32"},
    "eff_cAMP_vs_logRT":      {"r2": 0.74, "p": "<0.001"},
    "eff_cAMP_vs_logKi":      {"r2": 0.10, "p": "0.40"},
    "potency_vs_logRT":       {"r2": 0.077, "p": "0.44"},
    "kineticKD_vs_Ki":        {"r2": 0.99, "p": "<0.0001"},
}


def as_records():
    out = []
    for name, vals in COMPOUNDS.items():
        r = {"name": name}
        r.update(dict(zip(FIELDS, vals)))
        out.append(r)
    return out


def _check():
    """Internal consistency of the transcribed table.

    Guo et al. print koff rounded to two decimals but computed RT and the
    kinetic KD from unrounded rates, so 1/koff and koff/kon reproduce the
    printed RT and KD only to within rounding. The worst case in this table is
    LUF5550 at 7.4 %. A 12 % band therefore still catches any transcription
    error large enough to matter while tolerating the published rounding.
    """
    rs = as_records()
    assert len(rs) == 10
    worst_rt = worst_kd = 0.0
    for r in rs:
        rt = 1.0 / r["koff"]
        d_rt = abs(rt - r["RT"]) / r["RT"]
        assert d_rt < 0.12, (r["name"], "RT", rt, r["RT"])
        worst_rt = max(worst_rt, d_rt)

        kd = r["koff"] / r["kon"] * 1e9
        d_kd = abs(kd - r["kinetic_KD"]) / r["kinetic_KD"]
        assert d_kd < 0.20, (r["name"], "KD", kd, r["kinetic_KD"])
        worst_kd = max(worst_kd, d_kd)

    n_eff = sum(1 for r in rs if r["eff_impedance"] is not None)
    assert n_eff == 9, n_eff

    ids = [r["chembl_id"] for r in rs]
    assert len(set(ids)) == 10, "duplicate ChEMBL id"

    print(f"dataset self-check PASS: {len(rs)} compounds, {n_eff} with efficacy")
    print(f"  RT vs 1/koff        : max deviation {worst_rt*100:.1f} %")
    print(f"  kinetic KD vs koff/kon: max deviation {worst_kd*100:.1f} %")
    print("  all 10 ChEMBL ids distinct")


if __name__ == "__main__":
    _check()
