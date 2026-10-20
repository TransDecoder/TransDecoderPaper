# Genome annotation from assembled transcripts (GM12878)

Data and assemblies: `data/gm12878/README.md` (ENCODE ENCSR000AED rep1 short reads -> STAR + StringTie;
ENCSR838WFC rep2 PacBio FLNC -> minimap2 + LRAA; Trinity de novo -> minimap2). Reference: GENCODE v50,
MANE v1.5.

Scripts (`scripts/`):
- `run_td.sh <name> <fa> <structures> <abinitio|homology>`: PyTransDecoder LongOrfs (-S) + Predict
  (+ DIAMOND vs Swiss-Prot without human / Pfam for homology), then `genome-map` to genome CDS models.
- `eval_genome_cds.py`: CDS-structure evaluation against MANE Select and GENCODE complete CDSs
  (exact, compatible partial, assembled upper bound; prediction-centric precision).
- `start_labels.py`, `start_features.py`: label 5'-open ORFs by MANE (UTR extension vs truncated)
  with the start-rule features; used to compare keep-or-trim rule forms.
- `minlen_sweep.sh`: minimum ORF length sweep (90..50 aa).
