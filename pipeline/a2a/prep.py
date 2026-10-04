#!/usr/bin/env python
"""Prepare A2A receptors and the ten agonist ligands for docking.

Receptors: chain A protein atoms only. Co-crystal ligand, lipids and waters are
split out; the co-crystal ligand is kept separately as both the autobox
reference and the redocking validation target.

Ligands: SMILES from ChEMBL -> RDKit ETKDG conformer -> MMFF94 optimisation ->
Open Babel protonation at pH 7.4. No hand-drawn structures.
"""
import json, os, subprocess, sys, urllib.request, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dataset import COMPOUNDS, RECEPTORS

WORK = "/home/kjb9412/a2a_work"
STRUCT = f"{WORK}/structures"
LIG = f"{WORK}/ligands"
OBABEL = "/home/kjb9412/miniconda3/envs/dock/bin/obabel"


def fetch_smiles():
    """Pull canonical SMILES from ChEMBL for every compound in the set."""
    cache = f"{LIG}/smiles.json"
    if os.path.exists(cache):
        return json.load(open(cache))
    out = {}
    for name, vals in COMPOUNDS.items():
        cid = vals[0]
        url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{cid}.json"
        for _ in range(3):
            try:
                with urllib.request.urlopen(url, timeout=30) as r:
                    d = json.load(r)
                break
            except Exception:
                time.sleep(2)
        else:
            raise RuntimeError(f"ChEMBL fetch failed for {name} ({cid})")
        smi = (d.get("molecule_structures") or {}).get("canonical_smiles")
        if not smi:
            raise RuntimeError(f"no SMILES for {name} ({cid})")
        out[name] = {"chembl_id": cid, "smiles": smi,
                     "mw": (d.get("molecule_properties") or {}).get("full_mwt")}
        print(f"  {name:10s} {cid:15s} MW={out[name]['mw']}")
    os.makedirs(LIG, exist_ok=True)
    json.dump(out, open(cache, "w"), indent=1)
    return out


def split_receptor(pdb_id, lig_code):
    """Write protein-only receptor and the isolated co-crystal ligand."""
    src = f"{STRUCT}/{pdb_id.lower()}.pdb"
    rec_raw, lig_raw = f"{STRUCT}/{pdb_id}_rec_raw.pdb", f"{STRUCT}/{pdb_id}_native.pdb"
    with open(src) as f, open(rec_raw, "w") as fr, open(lig_raw, "w") as fl:
        for ln in f:
            if ln.startswith("ATOM") and ln[21] == "A":
                fr.write(ln)
            elif ln.startswith("HETATM") and ln[17:20].strip() == lig_code:
                fl.write(ln)
    fr_n = sum(1 for _ in open(rec_raw))
    fl_n = sum(1 for _ in open(lig_raw))
    if fl_n == 0:
        raise RuntimeError(f"{pdb_id}: no {lig_code} atoms found")

    rec = f"{STRUCT}/{pdb_id}_rec.pdb"
    subprocess.run([OBABEL, rec_raw, "-O", rec, "-h"], check=True,
                   capture_output=True)
    nat = f"{STRUCT}/{pdb_id}_native.sdf"
    subprocess.run([OBABEL, lig_raw, "-O", nat, "-h"], check=True,
                   capture_output=True)
    print(f"  {pdb_id}: {fr_n} protein atoms, {fl_n} {lig_code} atoms -> {rec}")
    return rec, nat


def build_ligand(name, smiles):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise RuntimeError(f"unparseable SMILES for {name}")
    m = Chem.AddHs(m)
    ps = AllChem.ETKDGv3()
    ps.randomSeed = 42
    if AllChem.EmbedMolecule(m, ps) != 0:
        raise RuntimeError(f"embedding failed for {name}")
    AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
    raw = f"{LIG}/{name}_raw.sdf"
    Chem.SDWriter(raw).write(m)
    out = f"{LIG}/{name}.sdf"
    subprocess.run([OBABEL, raw, "-O", out, "-p", "7.4"], check=True,
                   capture_output=True)
    chk = Chem.SDMolSupplier(out, removeHs=False)
    mol = next(iter(chk), None)
    if mol is None:
        raise RuntimeError(f"protonation produced unreadable sdf for {name}")
    fc = Chem.GetFormalCharge(mol)
    print(f"  {name:10s} atoms={mol.GetNumAtoms():3d} charge={fc:+d}")
    return out, fc


def main():
    os.makedirs(LIG, exist_ok=True)
    print("fetching SMILES from ChEMBL")
    smi = fetch_smiles()

    print("\nsplitting receptors")
    for pdb_id, meta in RECEPTORS.items():
        split_receptor(pdb_id, meta["ligand"])

    print("\nbuilding 3D ligands")
    charges = {}
    for name, rec in smi.items():
        _, fc = build_ligand(name, rec["smiles"])
        charges[name] = fc
    json.dump(charges, open(f"{LIG}/charges.json", "w"), indent=1)
    print(f"\nprepared {len(smi)} ligands, {len(RECEPTORS)} receptors")


if __name__ == "__main__":
    main()
