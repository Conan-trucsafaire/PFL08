# Phase 3a — resource-free coverage finding

**Date:** 2026-08-06
**Run:** `scripts/run_phase3a_offline.py`
**Input:** `data/frozen/02_SIGNED_EDGE_FIRST_PASS_CANDIDATES_FROZEN_2026-08-05.tsv` (53 records, SHA-256 `3E9E…52D6`)
**Resources consulted:** none. No effect value, z-score or p-value was read; no perturbation
dataset was downloaded. Outcome-blind by construction.

---

## Ruling

> **`EXPLORATORY_ONLY_ALREADY_DETERMINED`** — gate 1 cannot be met by the current frozen
> edge set under any resource coverage.

Protocol §3.3 requires **≥ 150 pooled T-visible testable edges**. Testable edges are a
subset of edges that survive the resource-independent exclusions, so those exclusions put a
hard ceiling on the gate that can be computed before a single byte is downloaded:

| ceiling | definition | count | vs 150 |
|---|---|---|---|
| **B1** | T-visible **and** structurally eligible, as the classifier stands | **24** | fails |
| **B2** | B1 + every `T_uncertain` edge promoted to T-visible by the §2.2 reviewer | **95** | **fails** |
| B3 | every edge except the intercellular ones | 213 | clears — but only by relaxing the detectability rule, which §2.3 prohibits |

**B2 is the number that settles it.** It already grants the human reviewer every visibility
promotion they could possibly make, and it is still 55 edges short. Reaching 150 would
additionally require assigning genes to pathway tokens (`mTORC1 activity`, `MAPK/ERK
activity`), to metabolites (`glycolysis`, `lactate`), or to non-coding RNAs not in the
screens' perturbation lists — each of which is "relaxing the detectability rule", named in
§2.3 as a prohibited remedy.

The remaining gates are not settled offline:

| gate | threshold | upper bound | verdict |
|---|---|---|---|
| `t_visible_testable_edges_pooled` | 150 | 95 | **CANNOT_PASS** |
| `completely_testable_cycles` | 25 | 26 | UNDECIDED |
| `contexts_with_50plus_testable_edges` | 2 | 2 | UNDECIDED |
| `nontransformed_contexts` | 1 | — | not boundable offline |
| `cycles_testable_in_2plus_contexts` | 10 | 26 | UNDECIDED |

Note how little headroom gate 2 has: **26 against a threshold of 25**, and that is an upper
bound assuming every provisionally-mapped symbol turns out to be real and present in both
primary contexts. Real coverage will only reduce it.

---

## Where the edges go

53 records → **70 circuits, 100% parsed** → **229 edges**, 70 cycles.

**Resource-independent attrition** (holds in every single-cell-line perturbation screen):

| reason | edges | share |
|---|---|---|
| structurally eligible | 110 | 48.0% |
| non-coding RNA endpoint | 73 | 31.9% |
| complex / pathway token endpoint | 25 | 10.9% |
| intercellular (TME, cell-type-prefixed) | 16 | 7.0% |
| metabolite or process endpoint | 5 | 2.2% |

**mRNA-visibility declaration** (protocol §3.2, visibility taken from the edge's *target*):

| stratum | edges | share |
|---|---|---|
| T_visible | 38 | 16.6% |
| T_uncertain | 162 | 70.7% |
| T_invisible | 29 | 12.7% |

The single largest sink is not the protein layer — it is **non-coding RNA endpoints, 31.9%
of all edges**. A third of this catalogue is miRNA / lncRNA / circRNA sponge circuits, and
those endpoints are neither perturbed nor measured in a CRISPRi Perturb-seq screen. That is
a structural property of the catalogue, not of the resources, and no alternative screen
fixes it.

**Cycle closure.** Only **39 of 70** circuits close at the symbol level — the first and last
node resolve to the same gene. The other 31 have an endpoint that is a non-coding RNA, a
TME node or a pathway token, so closure cannot even be stated, let alone tested. Of the 39,
**26** have every edge structurally eligible, and only **3** have every edge *both*
structurally eligible and T-visible (`39023169#0`, `39023169#1`, `41963301#0`). Estimand E2
— the closure statistic, the discriminative core of the whole design — would rest on those
3 cycles if the T-visible restriction were applied strictly.

---

## Two data-quality items for the §2.2 second review

1. **One curated positive loop has a negative sign product.** `42338478#0`:
   `STRAP ┤ DKK1-LRP6 inhibition → β-catenin/TCF4 → STRAP` multiplies to −1. The middle node
   is itself *named* as an inhibition, so the double negative is in the node name rather
   than in the arrow. A node-naming artefact, not a curation error — but every downstream
   loop-gain calculation reads the arrows, so it must be resolved before Phase 3b.
   Pinned by `tests/test_frozen_table.py::test_exactly_one_curated_loop_has_a_negative_sign_product`
   so a second instance cannot appear unnoticed.

2. **The 44.4% figure in protocol §3.2 is not reproduced by the classifier.** The protocol
   states that 44.4% of curated edges are non-transcriptional. The automated classifier
   yields 12.7% T_invisible and 70.7% T_uncertain — i.e. 83.4% *not affirmatively
   transcriptional*, and only 16.6% affirmatively so. The two figures are measuring
   different things and the discrepancy needs reconciling in the review, because §3.2's
   negative control (E8) depends on the T-invisible stratum being the *real* protein-layer
   set, not the residue of a conservative regex.

---

## What follows

The ruling is binding: **Phase 3b on the current frozen set produces no confirmatory claim.**
Protocol §2.3 permits exactly one remedy, and it is available:

> Extend signed-edge adjudication to the **248 un-adjudicated legacy records**, complete and
> re-freeze it before any effect value is read, then re-run Phase 3a.

Two things make this the right next move rather than a consolation prize:

- It is literature curation done **without reference to any perturbation outcome**, so it
  carries no leakage risk and needs no unmasking.
- The 53 adjudicated records yield 110 structurally eligible edges, ~2.1 per record. If the
  legacy records adjudicate at a similar rate and a similar visibility profile, roughly
  **90–120 further records** are needed to clear gate 1 at B2, or considerably more to clear
  it at B1. That is a scoping estimate from this run, not a promise — the legacy records
  have not been read, and their compartment and non-coding-RNA mix may differ.

Because the ceiling is set mostly by the catalogue's composition rather than by the screens,
one thing to decide **before** spending curation effort: whether the 248 legacy records are
richer in protein-coding, cell-intrinsic, transcriptional circuits than the 53 already
adjudicated. If they carry the same ~32% non-coding-RNA endpoint fraction, adjudicating all
248 still yields on the order of 100 additional B2-eligible edges — enough, but with no
margin. Sampling ~25 legacy records first and re-running this script would answer that for
a fraction of the effort.

**Independent of the gate**, and worth doing regardless of which branch is taken:

- `scripts/run_phase3a.py` still needs to run against R1/R2/R4/R6 to settle gates 2–5 and to
  produce the real (universe-verified) node resolution. The offline numbers are upper bounds.
- The §2.2 independent second review is still manuscript-blocking and still not done. The
  actionable artifact now exists: `results/phase3a_offline/phase3a_offline_NODE_REVIEW_WORKLIST.tsv`,
  **208 unique tokens** rather than 299 instances, with empty `reviewer_symbol`,
  `reviewer_visibility` and `reviewer_decision` columns to fill in.

## Honest framing of this result

This is a **coverage finding, not a result against the theory**. Nothing here says curated
positive feedback loops fail to close. It says the frozen 53-record set is too small and too
non-coding-heavy for genome-scale CRISPRi Perturb-seq to adjudicate closure at the
pre-registered evidentiary standard — which is precisely the kind of finding the gate was
put in front of the analysis to produce, and precisely the kind that gets buried when a
feasibility check is skipped.
