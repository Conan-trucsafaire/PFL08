# Governing documents and the clause → code map

The frozen protocol and charter live in Google Drive and are **not mirrored here**. A
frozen document copied into a second location can drift from the original, and a drifted
copy of a freeze is worse than a pointer to it. This file records their identity and maps
each operative clause onto the code that implements it.

## Governing documents

| Document | Location | Status |
|---|---|---|
| `10_ANALYSIS_PROTOCOL_PHASE3_PERTURBATION_RECIPROCITY_FROZEN_2026-08-05.md` | Drive `PFL08Claude版/`, file id `1gMhKZI5GESAwATLkmAO3ROaqUZoMU9Pn` | Frozen 2026-08-05, governing |
| `00_PROJECT_CHARTER_FROZEN_2026-08-05.md` | Drive `PFL生信08月/PFL_纯计算因果检验_2026-08-05/` | Frozen, still binding |
| `08_ANALYSIS_PROTOCOL_TOPOLOGY_STATISTICS_VALIDATION_FROZEN_2026-08-05.md` | same folder | **Not superseded** — withdrawal line continues under protocol §10 |
| `00_下一步战略_诊断与重新设计_2026-08-05.md` | Drive `PFL08Claude版/` | Strategy, four pillars |

## Frozen input

| File | SHA-256 |
|---|---|
| `data/frozen/02_SIGNED_EDGE_FIRST_PASS_CANDIDATES_FROZEN_2026-08-05.tsv` | `3E9E9CFB1E2922A2D94C5AAC83E7619199E377FD0E4C4932334DDCD5165452D6` |

Immutable (protocol §12.4). `tests/test_frozen_table.py` fails if it is edited.

## Clause → code

| Clause | Rule | Implemented in | Test |
|---|---|---|---|
| §2.1 | Direction and closure are never inferred from token order | `circuits.split_circuit` reads arrows only | `test_direction_is_never_inferred_from_token_order` |
| §2.2 | Independent second review is **manuscript-blocking** | worklist emitted; `second_human_review_completed: false` in every decision file | `test_every_record_is_still_pending_second_review` |
| §2.3 | Only permitted remedy for a failed gate is extending adjudication to the 248 legacy records | `gates.PERMITTED_REMEDY`, printed on every failure | — |
| §3 | Outcome-blind: no effect value read | `blind.BlindH5`, self-tested at startup; DepMap read with `nrows=0` | `test_matrix_groups_are_forbidden` |
| §3.1 | Testability: perturbed × measured × >500 kb apart × signed × mRNA-plausible | `testability.testable_in`, `CIS_WINDOW_BP` | `test_testable_requires_source_perturbed_and_target_measured` |
| §3.2 | mRNA-visibility declared before any lookup; visibility is a property of the **target** | `nodes.classify_visibility`, `graph.build_edges_and_cycles` | `test_edge_visibility_comes_from_the_target_not_the_source` |
| §3.2 | T-invisible edges are the prespecified **negative control**, not attrition | `testability.confirmatory_eligible` keeps them in the inventory | `test_t_invisible_edges_are_retained_as_the_negative_control` |
| §3.3 | Five gates, thresholds not editable | `gates.GATES` (immutable mapping) | `test_gate_thresholds_are_the_frozen_ones`, `test_gates_cannot_be_mutated` |
| §12.2 | Gate result written to a machine-readable lock before Phase 3b | `phase3a_gate_decision.json` | — |
| §12.4 | Raw inputs never modified | `data/frozen/` + checksum test | `test_frozen_table_is_unmodified` |
| §12.5 | Failed runs and null results are preserved, not deleted | `results/` is committed, including negative rulings | — |

## Not yet implemented here (Phase 3b, needs its own analysis lock)

§4 oriented edge score `w_c`, loop gain `g_c`; §5 comparators N1–N6; §6 estimands E1–E8;
§7 calibration and thresholds; §8.1 essential-gene stratification. None of it may be
written against real effect values before the §2.2 second review is complete and the
analysis lock is in place (§12.3 order of operations).

## Standing prohibitions (Charter §4, extended to abstract and title)

The words **necessary**, **sufficient**, **irreversible**, **prove**, **universal** — and
the phrases "proves positive feedback causes cancer", "hysteresis" — do not appear in any
claim, under any outcome.
