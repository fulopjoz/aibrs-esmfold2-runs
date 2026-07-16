# ESMFold2 runs for AibRS-Next (run on dichtator)

Local ESMFold2 on the 2x RTX 2080 Ti box is impractical (ESMC-6B is ~13 GB fp16, OOMs the 11 GB
Turing cards and does not auto-shard). This bundle runs on **dichtator** via the Singularity container
where ESMFold2 already works on GPU (~2 min/fold). Author the script here, run it there, push results back.

## Contents
- `targets.fasta` — the sequences to fold:
  - `tex1_mod2_adomain_1330-1903` (573 aa) and `tex1_mod13_adomain_13329-13882` (553 aa): the excised
    Tex1 (Q8NJX1) Aib-candidate A-domains. **B2 fold-check**: do they fold standalone (mean pLDDT, pTM)?
  - `alars_O58035_1-495_aminoacylation_module` (495 aa): the AlaRS aminoacylation module (ESMFold v1 gave
    1.64 A vs full-length; this re-checks it at ESMFold2 fidelity + reports the 7 pocket-residue pLDDTs).
- `fold_targets.py` — folds each target monomer, writes `results/<name>.cif` + `results/results.json`.
- `fold_with_ligand.py` — active-site co-fold: AlaRS module + ATP + 2 Mg + Aib (SMILES), vs Ala control.
  The ligand-aware experiment ESMFold2 enables. Honest caveat: a co-fold is a POSE model, not activation
  energetics (this project falsified binding as a ranker), so use it to visualise the pocket, not to rank.

## Run on dichtator
```bash
git clone https://github.com/fulopjoz/aibrs-esmfold2-runs.git
cd aibrs-esmfold2-runs
# fold-check (B2 + AlaRS module):
singularity exec --nv /home/hanzlv/share/esmfold2.sif python fold_targets.py
# optional active-site co-fold (Aib vs Ala):
singularity exec --nv /home/hanzlv/share/esmfold2.sif python fold_with_ligand.py
```

## Send results back (either works)
```bash
# git (preferred):
git add results/ && git commit -m "esmfold2 results" && git push
# or scp the results dir back to the analysis box:
scp -r results/ ubuntu@<box>:/home/ubuntu/ala_rs_aib/AlaRS-Aib-Validation-aibrs-next/esmfold2_dichtator/
```

## Notes
- ESMFold2 pLDDT is on a 0-1 scale; `fold_targets.py` reports it x100 (0-100).
- Defaults num_loops=8 / num_sampling_steps=100-150 (quality/speed balance; your example used 3/50).
- CCD ligands used: ATP, MG. Aib/Ala supplied as SMILES (Aib = `CC(C)(N)C(=O)O`).
