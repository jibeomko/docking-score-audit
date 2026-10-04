#!/usr/bin/env python
"""Do the docking scores predict affinity, kinetics, or efficacy?

Three questions, asked in order of increasing relevance to a real decision:

  Q1  score vs pKi         -- the quantity the scoring function is trained on
  Q2  score vs log RT      -- the kinetic quantity SPR measures
  Q3  score vs efficacy    -- the functional outcome

Guo et al. 2012 showed on this exact compound set that efficacy tracks log
residence time (r2 = 0.90) and not affinity (r2 = 0.13). Reproducing those two
numbers from the transcribed table is used here as a check that the table is
faithful, before any docking number is interpreted.
"""
import json, os, sys
import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dataset import as_records, PUBLISHED, RECEPTORS

WORK = os.environ.get("A2A_WORK", os.path.expanduser("~/a2a_work"))
OUT = f"{WORK}/out"


def corr(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 4:
        return None
    rho, p_s = stats.spearmanr(x, y)
    r, p_p = stats.pearsonr(x, y)
    return {"n": int(len(x)),
            "spearman": round(float(rho), 3), "spearman_p": round(float(p_s), 4),
            "pearson": round(float(r), 3), "pearson_p": round(float(p_p), 4),
            "r2": round(float(r) ** 2, 3)}


def boot_r2(x, y, n=5000, seed=0):
    """Percentile bootstrap CI on r2. n=10 needs the interval shown, not hidden."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 4:
        return None
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n):
        i = rng.integers(0, len(x), len(x))
        if len(set(x[i])) < 3:
            continue
        r = stats.pearsonr(x[i], y[i])[0]
        if np.isfinite(r):
            vals.append(r ** 2)
    if not vals:
        return None
    return [round(float(np.percentile(vals, 2.5)), 3),
            round(float(np.percentile(vals, 97.5)), 3)]


def main():
    recs = {r["name"]: r for r in as_records()}
    dock = json.load(open(f"{OUT}/dock_results.json"))

    # which receptors recovered their own co-crystal pose
    passed = {r["pdb"]: r.get("passed", False) for r in dock["redock"]}
    print("redocking validation")
    for r in dock["redock"]:
        print(f"  {r['pdb']} ({r['ligand']}, {r['resolution']} A): "
              f"rank1 RMSD {r.get('rmsd_rank1')} A -> "
              f"{'PASS' if r.get('passed') else 'FAIL'}")
    good = [p for p, v in passed.items() if v]
    print(f"  validated receptors: {good or 'none'}\n")

    # ---- per-compound score aggregation ----
    rows = []
    for name, rec in recs.items():
        row = {"name": name, "chembl_id": rec["chembl_id"],
               "Ki_nM": rec["Ki"], "pKi": -np.log10(rec["Ki"] * 1e-9),
               "koff": rec["koff"], "RT_min": rec["RT"],
               "logRT": np.log10(rec["RT"]),
               "kon": rec["kon"], "logkon": np.log10(rec["kon"]),
               "eff_impedance": rec["eff_impedance"],
               "eff_cAMP": rec["eff_cAMP"]}
        per_rec = {}
        for pdb in RECEPTORS:
            runs = [r for r in dock["prospective"][name].get(pdb, [])
                    if "cnn_affinity" in r]
            if not runs:
                continue
            cnn = [r["cnn_affinity"] for r in runs]
            vina = [r["vina"] for r in runs]
            per_rec[pdb] = {"cnn_mean": float(np.mean(cnn)),
                            "cnn_sd": float(np.std(cnn, ddof=1)) if len(cnn) > 1 else 0.0,
                            "vina_mean": float(np.mean(vina))}
            row[f"cnn_{pdb}"] = round(per_rec[pdb]["cnn_mean"], 3)
            row[f"cnn_sd_{pdb}"] = round(per_rec[pdb]["cnn_sd"], 3)
            row[f"vina_{pdb}"] = round(per_rec[pdb]["vina_mean"], 2)
        # ensemble summary: best score over validated receptors only
        pool = [v["cnn_mean"] for p, v in per_rec.items() if passed.get(p)] or \
               [v["cnn_mean"] for v in per_rec.values()]
        vpool = [v["vina_mean"] for p, v in per_rec.items() if passed.get(p)] or \
                [v["vina_mean"] for v in per_rec.values()]
        row["cnn_best"] = round(max(pool), 3) if pool else None
        row["cnn_mean_all"] = round(float(np.mean(pool)), 3) if pool else None
        row["vina_best"] = round(min(vpool), 2) if vpool else None
        row["seed_sd_max"] = round(max((v["cnn_sd"] for v in per_rec.values()),
                                       default=0.0), 3)
        rows.append(row)

    # ---- check the transcribed table against the published correlations ----
    eff = [r["eff_impedance"] for r in rows]
    lrt = [r["logRT"] for r in rows]
    lki = [np.log10(r["Ki_nM"]) for r in rows]
    checks = {
        "eff_impedance_vs_logRT": corr(lrt, eff),
        "eff_impedance_vs_logKi": corr(lki, eff),
        "eff_cAMP_vs_logRT": corr(lrt, [r["eff_cAMP"] for r in rows]),
        "eff_cAMP_vs_logKi": corr(lki, [r["eff_cAMP"] for r in rows]),
    }
    print("table fidelity check -- reproducing Guo et al. 2012 from the transcribed values")
    for k, v in checks.items():
        pub = PUBLISHED.get(k, {}).get("r2")
        if v:
            flag = "ok" if pub is None or abs(v["r2"] - pub) <= 0.15 else "DRIFT"
            print(f"  {k:26s} r2={v['r2']:.2f}  published={pub}  [{flag}]")

    # ---- the actual question ----
    targets = {"pKi": [r["pKi"] for r in rows],
               "logRT": lrt,
               "logkoff": [np.log10(r["koff"]) for r in rows],
               "logkon": [r["logkon"] for r in rows],
               "eff_impedance": eff,
               "eff_cAMP": [r["eff_cAMP"] for r in rows]}
    preds = {"CNNaffinity_best": [r["cnn_best"] for r in rows],
             "CNNaffinity_mean": [r["cnn_mean_all"] for r in rows],
             "Vina_best": [r["vina_best"] for r in rows]}

    table = {}
    print("\ndocking score vs experiment")
    for pn, pv in preds.items():
        table[pn] = {}
        for tn, tv in targets.items():
            c = corr(pv, tv)
            if c:
                c["r2_ci95"] = boot_r2(pv, tv)
                table[pn][tn] = c
                print(f"  {pn:18s} vs {tn:14s} "
                      f"rho={c['spearman']:+.3f} (p={c['spearman_p']:.3f})  "
                      f"r2={c['r2']:.3f} CI{c['r2_ci95']}")

    json.dump({"rows": rows, "fidelity_checks": checks,
               "published": PUBLISHED, "score_vs_experiment": table,
               "redock": dock["redock"], "meta": dock["meta"],
               "sampling_sweep": dock.get("sampling_sweep", [])},
              open(f"{OUT}/analysis.json", "w"), indent=1, default=float)

    import csv
    keys = sorted({k for r in rows for k in r})
    order = ["name", "chembl_id", "Ki_nM", "pKi", "koff", "RT_min", "logRT",
             "kon", "eff_impedance", "eff_cAMP", "cnn_best", "cnn_mean_all",
             "vina_best", "seed_sd_max"]
    order += [k for k in keys if k not in order]
    with open(f"{OUT}/results.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=order)
        w.writeheader()
        for r in sorted(rows, key=lambda r: -r["RT_min"]):
            w.writerow(r)
    print(f"\nwrote {OUT}/analysis.json and {OUT}/results.csv")


def selection_simulation():
    """If you picked compounds by score, what would you have picked?

    The correlation coefficients above are the statistician's view. This is the
    medicinal chemist's view: rank the set by each available criterion, take the
    top three, and report the mean functional efficacy of that selection. It
    converts a weak correlation into the decision it would actually drive.
    """
    import json
    import numpy as np
    A = json.load(open(f"{OUT}/analysis.json"))
    rows = [r for r in A["rows"] if r["eff_impedance"] is not None]
    k = 3
    criteria = {
        "residence time (SPR)":   lambda r: r["RT_min"],
        "affinity (1/Ki)":        lambda r: -r["Ki_nM"],
        "GNINA CNNaffinity":      lambda r: r["cnn_best"],
        "Vina (most negative)":   lambda r: -r["vina_best"],
    }
    base = float(np.mean([r["eff_impedance"] for r in rows]))
    best = float(np.mean(sorted((r["eff_impedance"] for r in rows),
                                reverse=True)[:k]))
    out = {"n": len(rows), "top_k": k,
           "mean_efficacy_whole_set": round(base, 1),
           "mean_efficacy_perfect_pick": round(best, 1), "by_criterion": {}}
    print(f"\nselection simulation: top {k} of {len(rows)} by each criterion")
    print(f"  {'criterion':24s} {'picks':34s} mean efficacy")
    for label, key in criteria.items():
        picks = sorted(rows, key=key, reverse=True)[:k]
        m = float(np.mean([p["eff_impedance"] for p in picks]))
        names = ", ".join(p["name"] for p in picks)
        out["by_criterion"][label] = {"picks": [p["name"] for p in picks],
                                      "mean_efficacy": round(m, 1),
                                      "vs_random": round(m - base, 1)}
        print(f"  {label:24s} {names:34s} {m:5.1f}  ({m - base:+.1f} vs random)")
    print(f"  {'(random pick)':24s} {'-':34s} {base:5.1f}")
    print(f"  {'(perfect pick)':24s} {'-':34s} {best:5.1f}")
    A["selection_simulation"] = out
    json.dump(A, open(f"{OUT}/analysis.json", "w"), indent=1, default=float)
    return out


if __name__ == "__main__":
    main()
    selection_simulation()
