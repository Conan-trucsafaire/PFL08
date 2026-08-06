# Amendment A1 — node-resolution defect in the first-pass classifier

**Date:** 2026-08-06
**Made before or after which outcomes were read:** **Before any outcome.** No perturbation
resource was consulted; no effect value, z-score or p-value exists in this project yet.
The amendment is therefore outcome-blind in the sense of protocol §2.3 and §12.3.
**Nature:** correction of a matching defect. No threshold, estimand, comparator or
curation vocabulary was changed.

## What was wrong

The first-pass classifier (`11_phase3a.py`, Drive `PFL08Claude版/`) decided whether a node
token named a pathway or complex with a **bare prefix test**:

```python
for c in COMPLEX_TOKENS:
    if t.lower().startswith(c.lower()) or f" {c.lower()} " in f" {t.lower()} ":
        # -> complex_or_pathway, no gene assigned
```

`COMPLEX_TOKENS` contains `"Wnt"`, `"AKT"`, `"MAPK"`, `"ERK"`, `"PI3K"`, `"IKK"`, `"p38"`.
Any real gene whose symbol *begins with* one of those strings was therefore swallowed:

| token | classified as | should be |
|---|---|---|
| `WNT7B` | `complex_or_pathway` | gene `WNT7B` |
| `secreted WNT7B` | gene `WNT7B` (no prefix — inconsistent with the line above) | gene `WNT7B` |
| `AKT1` | `complex_or_pathway` | gene `AKT1` |
| `MAPK14` | `complex_or_pathway` | gene `MAPK14` |

Two consequences, both in the direction of under-counting:

1. Those edges were reported as permanently untestable, when the gene is perturbed in
   Replogle K562 and measured in every context.
2. **Closure detection broke.** `secreted WNT7B -| FTO activity -| TCF7L2 abundance/activity -> WNT7B`
   was read as *not closed*, because the bare `WNT7B` endpoint had become a pathway with
   no symbol while the `secreted WNT7B` endpoint had not.

A second, smaller defect: `SYMBOL_RE`'s optional `[A-Z0-9-]{0,4}` tail is greedy across
hyphens and `\b` is satisfied immediately after one, so `IL-1β` yielded the "symbol" `IL-`
and `DKK1-LRP6 inhibition` yielded `DKK1-`. Neither is a gene symbol, and each occupied
the node's symbol slot, blocking the correct one.

## What was changed

`pfl_phase3/nodes.py` only:

- complex/pathway tokens are matched with a **token boundary on both sides**
  (`(?<![A-Za-z0-9])TOKEN(?![A-Za-z0-9])`), so `Wnt` still matches `Wnt signaling` and
  `AKT` still matches `AKT activity`, but neither matches `WNT7B` or `AKT1`;
- symbol hits are stripped of trailing hyphens and dropped below two characters.

**The `COMPLEX_TOKENS`, `PROCESS_TOKENS`, `T_VISIBLE_KEYS` and `T_INVISIBLE_KEYS`
vocabularies are the author's curation and were not touched.** Only the matching
mechanics changed. Regression tests: `tests/test_nodes.py::test_gene_starting_with_a_pathway_name_is_not_a_pathway`,
`::test_pathway_matching_still_fires_on_whole_tokens`, `::test_trailing_hyphen_is_not_part_of_a_symbol`.

## Effect on the Phase 3a numbers

History is not overwritten (§12.4). Both runs, on the same frozen input:

| quantity | before A1 | after A1 |
|---|---|---|
| circuits parsed | 70 / 70 | 70 / 70 |
| edges | 229 | 229 |
| nodes classified `complex_or_pathway` | 16 | 15 |
| nodes classified `single_symbol` | 177 | 179 |
| cycles CLOSED | 37 | **39** |
| edges structurally eligible | 110 | 110 |
| gate 1 ceiling B1 (T-visible & eligible) | 24 | 24 |
| gate 2 ceiling (closed & fully eligible cycles) | 25 | **26** |

**The ruling does not change.** Gate 1 remains `CANNOT_PASS` before and after, by a wide
margin, so the defect was not what determined the Phase 3a outcome. It is recorded because
it would have silently corrupted Phase 3b edge scoring, and because the closure count —
the input to estimand E2 — moved.

## Still open

`MAPK14/p38 activity` remains `complex_or_pathway`, because `p38` is in the curated
pathway vocabulary. MAPK14 *is* p38α, a real perturbable gene. Whether to promote this
node is a curation judgement for the §2.2 second reviewer, not a code change, and it is
flagged in `results/phase3a_offline/phase3a_offline_MANUAL_REVIEW_WORKLIST.tsv`.
