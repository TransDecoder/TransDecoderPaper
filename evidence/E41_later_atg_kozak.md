# E41. Later in-frame ATGs scored by start context, release code (2026-11-02)

Replaces the development-code comparison (PyTD efc1450 `--start-atg-choice pwm`, `pytd_dev_pwm`) that supported
the Methods sentence on considering only the first in-frame ATG, so the claim rests on the release code
(Brian Haas, 2026-11-02, who proposed scoring the alternative starts directly). The "PWM, scale corrected"
column (PyTD 27bfb5e) and its Results sentence (median exact F1 0.728) were dropped at the same time (BH: not
needed).

## Claims supported
Methods, "Start-codon refinement" (sentence after "Only the first in-frame ATG is considered");
Supplementary Table S6b. Consistent with Results "Remaining start errors" (C. reinhardtii later starts).

## Method
`analyses/start_param_sensitivity/later_atg_kozak_v1.py` -> `data/start_refinement/later_atg_kozak.tsv`.
- Inputs: the stranded `complete` input of each organism; complete-type truth CDSs with an ATG start.
- Model: `start_refinement.kozak.tsv` of the pinned release run on that input
  (`$P/work/release_pin/bench/<org>/complete/rel/`, PyTD d3892b5), scored with `pytransdecoder.core.kozak.score`
  (bits), i.e. the start-context model the release trains and uses.
- Candidates: in-frame ATGs from the first one after the upstream in-frame stop (or transcript start) up to the
  annotated stop (stop set TGA only for T. thermophila).
- Rule: move from the first ATG to the best-scoring later ATG within W codons of it (10, 50, all) when it scores
  at least M bits above the first ATG (M = 0-4; M = 0 means strictly higher). gain = annotated later start chosen;
  loss = annotated first-ATG start moved away from; later_wrong = moved to a different later ATG.
- Runtime: about 2 min, one core.

## Results
- Annotated start at a later ATG: 1.0-4.9% of CDSs in seven organisms, 19.5% in C. reinhardtii (3,438 of 17,621).
- Of CDSs starting at the first ATG, a later ATG within 50 codons scores higher in 10.5-20.9% (within 10 codons
  2.6-7.7%).
- Windows 10 or 50 codons, margins 2 or 4 bits: net loss of exact matches at every setting in all seven other
  organisms (e.g. human -344, -951, -205 at 10/2, 50/2, 50/4). C. reinhardtii: +58, +20, +83; best 50 codons and
  4 bits, 247 gained, 164 lost (7% of its later-ATG starts recovered). Margin 0 or no window: net loss everywhere,
  including C. reinhardtii.
- Agrees with the development-code result it replaces (PWM choice +69 in C. reinhardtii, lower in the others
  except P. falciparum +1; E09).

## Checks
Notebook 03, section "Later ATGs scored by start context" (`check_text` rows). The "later ATG" fractions
match `data/start_refinement/later_atg.tsv` for C. reinhardtii (19.6% among transcripts with annotated 5' UTRs
there; 19.5% of all complete-type CDSs here, different denominator).

## Caveats
- An upper bound on what start context alone can do with this model; a rule combining context with distance or
  coding score was not tested.
- The Kozak model is trained on the input's own longest ORFs, which include many of these CDSs, as in the release.
