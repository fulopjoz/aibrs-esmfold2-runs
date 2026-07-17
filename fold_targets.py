#!/usr/bin/env python
"""ESMFold2 fold-check of the AibRS-Next targets (run on dichtator via the esmfold2.sif container).

Folds each sequence in targets.fasta as a standalone monomer, writes <id>.cif, and records mean pLDDT
+ pTM + (for the AlaRS module) the 7 pocket-residue pLDDTs into results.json.

Run:  singularity exec --nv /home/hanzlv/share/esmfold2.sif python fold_targets.py
(ESMFold2 pLDDT is on a 0-1 scale; this script reports it x100.)
"""
import json, time
from pathlib import Path
from esm.models.esmfold2 import ESMFold2InputBuilder, ProteinInput, StructurePredictionInput
from transformers.models.esmfold2.modeling_esmfold2 import ESMFold2Model

HERE = Path(__file__).resolve().parent
OUT = HERE / "results"; OUT.mkdir(exist_ok=True)
# Per-target pocket positions (1-based); report per-residue pLDDT for these when the name matches.
POCKETS = {
    "1-495": [192, 193, 213, 215, 216, 217, 249],           # AlaRS O58035 pocket (no offset)
    "chutrakul": [211, 212, 215, 252, 273, 275, 296, 304],  # Chutrakul A-domain Stachelhaus code (excised coords)
}

def read_fasta(path):
    recs, cur, seq = [], None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if cur: recs.append((cur, "".join(seq)))
            cur = line[1:].strip(); seq = []
        elif line.strip():
            seq.append(line.strip())
    if cur: recs.append((cur, "".join(seq)))
    return recs

def main():
    print("[esmfold2] loading biohub/ESMFold2 (GPU)...", flush=True)
    t0 = time.time()
    model = ESMFold2Model.from_pretrained("biohub/ESMFold2").cuda().eval()
    print(f"[esmfold2] loaded in {time.time()-t0:.0f}s", flush=True)
    builder = ESMFold2InputBuilder()
    results = []
    for name, seq in read_fasta(HERE / "targets.fasta"):
        print(f"[fold] {name} ({len(seq)} aa)...", flush=True)
        tf = time.time()
        spi = StructurePredictionInput(sequences=[ProteinInput(id="A", sequence=seq)])
        res = builder.fold(model, spi, num_loops=8, num_sampling_steps=100,
                           num_diffusion_samples=1, seed=0)
        (OUT / f"{name}.cif").write_text(res.complex.to_mmcif())
        plddt = res.plddt.detach().float().cpu().numpy().ravel() * 100.0  # 0-1 -> 0-100
        rec = {"name": name, "length": len(seq), "seconds": round(time.time()-tf, 1),
               "mean_plddt": round(float(plddt.mean()), 1),
               "frac_plddt_ge70": round(float((plddt >= 70).mean()), 3),
               "ptm": round(float(res.ptm), 3), "cif": f"results/{name}.cif"}
        # per-residue pLDDT at the pocket (AlaRS module or Chutrakul A-domain)
        for _key, _positions in POCKETS.items():
            if _key in name and plddt.size >= max(_positions):
                rec["pocket_plddt"] = {p: round(float(plddt[p-1]), 1) for p in _positions}
                break
        results.append(rec)
        print(f"[fold]   mean pLDDT={rec['mean_plddt']} pTM={rec['ptm']} ({rec['seconds']}s)", flush=True)
    (OUT / "results.json").write_text(json.dumps(results, indent=2))
    print(f"[esmfold2] DONE -> {OUT}/results.json + per-target .cif", flush=True)
    for r in results:
        print(f"  {r['name']:45} L={r['length']:4} meanpLDDT={r['mean_plddt']:5} pTM={r['ptm']}")

if __name__ == "__main__":
    main()
