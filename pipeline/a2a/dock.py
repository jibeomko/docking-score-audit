#!/usr/bin/env python
"""Dock the ten A2A agonists into three agonist-bound A2A structures.

Two stages, in this order, so that validation cannot be tuned after the fact:

1. Redocking validation. Each co-crystal ligand is docked back into its own
   structure and the pose RMSD against the crystal coordinates is recorded.
   A structure whose own ligand is not recovered is flagged in the output and
   the analysis stage down-weights it.

2. Prospective docking. All ten agonists into all three receptors, three seeds
   each, so that score spread across seeds is measurable rather than assumed.
"""
import json, os, subprocess, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dataset import COMPOUNDS, RECEPTORS

WORK = os.environ.get("A2A_WORK", os.path.expanduser("~/a2a_work"))
GNINA = os.environ.get("GNINA", "gnina")
CUDNN = os.environ.get("CUDNN_LIB", "")   # directory with libcudnn, if gnina needs it
SEEDS = (42, 7, 1234)
EXHAUSTIVENESS = 16
NUM_MODES = 20
SWEEP = (16, 32, 64)   # sampling-adequacy check, stage 0


def env():
    e = dict(os.environ)
    if CUDNN:
        e["LD_LIBRARY_PATH"] = CUDNN + ":" + e.get("LD_LIBRARY_PATH", "")
    return e


def run_gnina(rec, lig, box_ref, out, seed):
    cmd = [GNINA, "-r", rec, "-l", lig,
           "--autobox_ligand", box_ref, "--autobox_add", "6",
           "-o", out, "--seed", str(seed),
           "--num_modes", str(NUM_MODES),
           "--exhaustiveness", str(EXHAUSTIVENESS)]
    t0 = time.time()
    p = subprocess.run(cmd, env=env(), capture_output=True, text=True)
    if p.returncode != 0 or not os.path.exists(out):
        return None, p.stderr[-400:], time.time() - t0
    return parse_sdf(out), None, time.time() - t0


def parse_sdf(path):
    """Read CNNaffinity / CNNscore / Vina affinity for every pose."""
    from rdkit import Chem
    poses = []
    for m in Chem.SDMolSupplier(path, removeHs=False, sanitize=False):
        if m is None:
            continue
        g = lambda k: float(m.GetProp(k)) if m.HasProp(k) else None
        poses.append({"cnn_affinity": g("CNNaffinity"),
                      "cnn_score": g("CNNscore"),
                      "vina": g("minimizedAffinity")})
    return poses


def best_rms(ref_sdf, probe_sdf, pose_idx=0):
    """Heavy-atom RMSD of one docked pose against the crystal coordinates.

    GNINA writes poses of the input molecule with atom order preserved, and the
    redocking input here *is* the crystal ligand, so a direct positional RMSD
    over heavy atoms is the correct comparison and needs no atom mapping.
    Ligands lifted straight out of a PDB file often carry valences RDKit will
    not sanitise, so everything is read with sanitize=False. A symmetry-aware
    RMSD is attempted first and used when the molecule does sanitise, since it
    is the fairer number for rings and equivalent substituents.
    """
    import numpy as np
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign

    ref = next(iter(Chem.SDMolSupplier(ref_sdf, removeHs=False, sanitize=False)), None)
    probes = [m for m in Chem.SDMolSupplier(probe_sdf, removeHs=False,
                                            sanitize=False) if m is not None]
    if ref is None or pose_idx >= len(probes):
        return None
    probe = probes[pose_idx]

    # preferred: symmetry-aware, only if both molecules sanitise cleanly
    try:
        r = Chem.Mol(ref); p = Chem.Mol(probe)
        Chem.SanitizeMol(r); Chem.SanitizeMol(p)
        r = Chem.RemoveHs(r); p = Chem.RemoveHs(p)
        return round(float(rdMolAlign.CalcRMS(p, r)), 3)
    except Exception:
        pass

    # fallback: positional RMSD over heavy atoms, atom order assumed preserved
    hr = [a.GetIdx() for a in ref.GetAtoms() if a.GetAtomicNum() > 1]
    hp = [a.GetIdx() for a in probe.GetAtoms() if a.GetAtomicNum() > 1]
    if len(hr) != len(hp) or not hr:
        return None
    a = ref.GetConformer().GetPositions()[hr]
    b = probe.GetConformer().GetPositions()[hp]
    return round(float(np.sqrt(((a - b) ** 2).sum(axis=1).mean())), 3)


def sampling_sweep(S):
    """Is a failed redocking a search failure or a scoring failure?

    Each co-crystal ligand is redocked at increasing exhaustiveness. If the
    best recoverable pose stops improving while the budget keeps growing, the
    crystal pose is not the scoring optimum and no amount of extra search will
    find it. That distinction decides whether the prospective runs below are
    adequately sampled at exhaustiveness 16.
    """
    rows = []
    for pdb, meta in RECEPTORS.items():
        for ex in SWEEP:
            out = f"{WORK}/out/sweep_{pdb}_ex{ex}.sdf"
            cmd = [GNINA, "-r", f"{S}/{pdb}_rec.pdb", "-l", f"{S}/{pdb}_native.sdf",
                   "--autobox_ligand", f"{S}/{pdb}_native.pdb", "--autobox_add", "6",
                   "-o", out, "--seed", "42", "--num_modes", str(NUM_MODES),
                   "--exhaustiveness", str(ex)]
            t0 = time.time()
            r = subprocess.run(cmd, env=env(), capture_output=True, text=True)
            dt = time.time() - t0
            if r.returncode != 0 or not os.path.exists(out):
                rows.append({"pdb": pdb, "exhaustiveness": ex, "error": True})
                print(f"  {pdb} ex={ex}: FAILED", flush=True)
                continue
            poses = parse_sdf(out)
            rs = [best_rms(f"{S}/{pdb}_native.sdf", out, i)
                  for i in range(min(len(poses), NUM_MODES))]
            rs = [x for x in rs if x is not None]
            best = min(rs) if rs else None
            rows.append({"pdb": pdb, "exhaustiveness": ex,
                         "rmsd_rank1": rs[0] if rs else None,
                         "rmsd_best": best,
                         "rank_of_best": (rs.index(best) + 1) if best else None,
                         "cnn_affinity_rank1": poses[0]["cnn_affinity"],
                         "seconds": round(dt, 1)})
            print(f"  {pdb} ex={ex:3d}: rank1={rs[0] if rs else None} A  "
                  f"best={best} A  CNNaff={poses[0]['cnn_affinity']:.3f}  "
                  f"{dt:.0f}s", flush=True)
    return rows


def main():
    os.makedirs(f"{WORK}/out", exist_ok=True)
    S, L = f"{WORK}/structures", f"{WORK}/ligands"
    OUT_PARTIAL = f"{WORK}/out/results_partial.json"
    results = {"meta": {"seeds": list(SEEDS), "exhaustiveness": EXHAUSTIVENESS,
                        "num_modes": NUM_MODES, "gnina": "v1.3.3",
                        "autobox_add": 6,
                        "box_reference": "co-crystal ligand of each structure"},
               "redock": [], "prospective": {}}

    print("=== stage 0: sampling adequacy sweep ===", flush=True)
    results["sampling_sweep"] = sampling_sweep(S)
    json.dump(results, open(f"{OUT_PARTIAL}", "w"), indent=1)

    print("\n=== stage 1: redocking validation ===", flush=True)
    for pdb, meta in RECEPTORS.items():
        rec, nat = f"{S}/{pdb}_rec.pdb", f"{S}/{pdb}_native.sdf"
        out = f"{WORK}/out/redock_{pdb}.sdf"
        poses, err, dt = run_gnina(rec, nat, f"{S}/{pdb}_native.pdb", out, 42)
        if poses is None:
            print(f"  {pdb}: FAILED {err}", flush=True)
            results["redock"].append({"pdb": pdb, "error": err})
            continue
        rms1 = best_rms(nat, out, 0)
        rms_all = [best_rms(nat, out, i) for i in range(min(len(poses), NUM_MODES))]
        rms_all = [x for x in rms_all if x is not None]
        best = min(rms_all) if rms_all else None
        rank = (rms_all.index(best) + 1) if best is not None else None
        results["redock"].append({
            "pdb": pdb, "ligand": meta["ligand"], "native": meta["native"],
            "resolution": meta["resolution"], "n_poses": len(poses),
            "rmsd_rank1": rms1, "rmsd_best": best, "rank_of_best": rank,
            "cnn_affinity_rank1": poses[0]["cnn_affinity"],
            "vina_rank1": poses[0]["vina"],
            "passed": (rms1 is not None and rms1 <= 2.0),
            "seconds": round(dt, 1)})
        print(f"  {pdb} ({meta['ligand']}): rank1 RMSD={rms1} best={best} "
              f"@rank{rank}  CNNaff={poses[0]['cnn_affinity']}  {dt:.0f}s", flush=True)

    json.dump(results, open(f"{WORK}/out/results_partial.json", "w"), indent=1)

    print("\n=== stage 2: prospective docking ===", flush=True)
    total = len(COMPOUNDS) * len(RECEPTORS) * len(SEEDS)
    done = 0
    for name in COMPOUNDS:
        results["prospective"][name] = {}
        for pdb, meta in RECEPTORS.items():
            runs = []
            for seed in SEEDS:
                out = f"{WORK}/out/{name}__{pdb}__s{seed}.sdf"
                poses, err, dt = run_gnina(
                    f"{S}/{pdb}_rec.pdb", f"{L}/{name}.sdf",
                    f"{S}/{pdb}_native.pdb", out, seed)
                done += 1
                if poses is None:
                    print(f"  [{done}/{total}] {name}/{pdb}/s{seed} FAILED", flush=True)
                    runs.append({"seed": seed, "error": err})
                    continue
                runs.append({"seed": seed,
                             "cnn_affinity": poses[0]["cnn_affinity"],
                             "cnn_score": poses[0]["cnn_score"],
                             "vina": poses[0]["vina"],
                             "n_poses": len(poses), "seconds": round(dt, 1)})
                print(f"  [{done}/{total}] {name:10s} {pdb} s{seed:<5d} "
                      f"CNNaff={poses[0]['cnn_affinity']:.3f} "
                      f"vina={poses[0]['vina']:.2f}  {dt:.0f}s", flush=True)
            results["prospective"][name][pdb] = runs
            json.dump(results, open(f"{WORK}/out/results_partial.json", "w"), indent=1)

    json.dump(results, open(f"{WORK}/out/dock_results.json", "w"), indent=1)
    print("\nwrote out/dock_results.json", flush=True)


if __name__ == "__main__":
    main()
