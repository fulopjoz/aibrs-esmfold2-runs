#!/usr/bin/env python
"""ESMFold2 active-site co-fold: AlaRS aminoacylation module + ATP + Aib + Mg (the pre-adenylation state).

This is the ligand-aware experiment ESMFold2 enables (SMILES + CCD ligands). It co-folds the AlaRS
1-495 module with ATP (CCD), catalytic Mg (CCD), and Aib (SMILES) to inspect whether Aib is
accommodated in a near-attack geometry to the ATP alpha-phosphate. HONEST CAVEAT: a co-fold is a
plausibility/pose model, NOT activation energetics; binding pose does not equal activation (this
project falsified binding as a ranker). Use it to visualise the pocket, not to rank.

Run:  singularity exec --nv /home/hanzlv/share/esmfold2.sif python fold_with_ligand.py
"""
import json, time
from pathlib import Path
from esm.models.esmfold2 import (ESMFold2InputBuilder, ProteinInput, LigandInput,
                                 StructurePredictionInput)
from transformers.models.esmfold2.modeling_esmfold2 import ESMFold2Model

HERE = Path(__file__).resolve().parent
OUT = HERE / "results"; OUT.mkdir(exist_ok=True)
AIB_SMILES = "CC(C)(N)C(=O)O"          # 2-aminoisobutyric acid
ALA_SMILES = "C[C@@H](N)C(=O)O"        # L-alanine control

def _iter():
    cur, seq = None, []
    for line in (HERE / "targets.fasta").read_text().splitlines():
        if line.startswith(">"):
            if cur: yield cur, seq
            cur, seq = line[1:].strip(), []
        elif line.strip(): seq.append(line.strip())
    if cur: yield cur, seq

def main():
    seq = next(("".join(s) for h, s in _iter() if "1-495" in h), None)
    if not seq:
        raise SystemExit("AlaRS 1-495 target not found in targets.fasta")
    print("[esmfold2] loading model (GPU)...", flush=True)
    model = ESMFold2Model.from_pretrained("biohub/ESMFold2").cuda().eval()
    builder = ESMFold2InputBuilder()
    out = {}
    for tag, aa_smiles in (("aib", AIB_SMILES), ("ala", ALA_SMILES)):
        print(f"[cofold] AlaRS module + ATP + Mg + {tag} ...", flush=True)
        t = time.time()
        spi = StructurePredictionInput(sequences=[
            ProteinInput(id="A", sequence=seq),
            LigandInput(id="ATP", ccd=["ATP"]),
            LigandInput(id="MG1", ccd=["MG"]),
            LigandInput(id="MG2", ccd=["MG"]),
            LigandInput(id="AA", smiles=aa_smiles),
        ])
        res = builder.fold(model, spi, num_loops=8, num_sampling_steps=150,
                           num_diffusion_samples=1, seed=0)
        (OUT / f"alars_active_site_{tag}.cif").write_text(res.complex.to_mmcif())
        out[tag] = {"mean_plddt": round(float(res.plddt.mean()) * 100, 1),
                    "ptm": round(float(res.ptm), 3), "iptm": round(float(res.iptm), 3),
                    "seconds": round(time.time() - t, 1),
                    "cif": f"results/alars_active_site_{tag}.cif"}
        print(f"[cofold]   {tag}: pLDDT={out[tag]['mean_plddt']} ipTM={out[tag]['iptm']} "
              f"({out[tag]['seconds']}s)", flush=True)
    (OUT / "active_site_cofold.json").write_text(json.dumps(out, indent=2))
    print("[esmfold2] DONE -> results/active_site_cofold.json + .cif (inspect the pocket in PyMOL)", flush=True)

if __name__ == "__main__":
    main()
