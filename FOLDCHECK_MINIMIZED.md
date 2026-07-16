# Fold-check the conservation-minimized AlaRS constructs

From `results/v3/aminoacylation_redesign/domain_minimization.py` (717-AlaRS MSA conservation):
- `alars_module_1-495` (495 aa) — patent construct, control (ESMFold v1/2 ~94 pLDDT).
- `alars_module_minus_insertions_384aa` — 1-495 minus the 111 aa of O58035-specific insertions the
  family lacks (N-terminal 1-57 + internal loops 154-159/223-234/285-296/310-336). Lowest-risk minimization.
- `alars_conserved_core_21-489` (469 aa) — conserved catalytic core.

Run on dichtator to test whether the minimized construct still folds (junctions from deleted insertions
are the risk). Edit fold_targets.py to read minimized_constructs.fasta, or:
  singularity exec --nv /home/hanzlv/share/esmfold2.sif python fold_targets.py
(after pointing its `read_fasta` at minimized_constructs.fasta). Report mean pLDDT + pTM + pocket pLDDT;
a big drop at a deleted-insertion junction means that insertion is structurally load-bearing (keep it).
